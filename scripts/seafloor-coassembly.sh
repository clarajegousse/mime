#!/usr/bin/env bash
# seafloor-coassembly

conda activate anvio-6.2

# Set you own project directory
MIME_DIR=/users/work/cat3/projects/mime

# select clean reads (sunbeam output)
R1s=$(ls /users/work/cat3/projects/mime/sunbeam_output/qc/decontam/*_wgs_1.fastq.gz | grep "wgs" | grep -e 1045m -e 1004m -e 90m -e 470m -e MC | python -c 'import sys; print( ",".join([x.strip() for x in sys.stdin.readlines()]))' )
R2s=$(ls /users/work/cat3/projects/mime/sunbeam_output/qc/decontam/*_wgs_2.fastq.gz | grep "wgs" | grep -e 1045m -e 1004m -e 90m -e 470m -e MC | python -c 'import sys; print( ",".join([x.strip() for x in sys.stdin.readlines()]))' )

ASB=coassembly_wgs_bottom

MIN_CONTIG_SIZE=1000
NUM_THREADS=12
ASSEMBLY_DIR=$MIME_DIR/results/coassembly/$ASB
PROFILE_DIR=$MIME_DIR/results/coassembly/$ASB/profile
MERGE_DIR=$MIME_DIR/results/coassembly/$ASB/merge
SUMMARY_DIR=$MIME_DIR/results/coassembly/$ASB/summary
  
megahit -1 $R1s -2 $R2s --min-contig-len $MIN_CONTIG_SIZE -m 0.85 -o $ASSEMBLY_DIR/ -t $NUM_THREADS

cd $ASSEMBLY_DIR

# ---------

anvi-script-reformat-fasta $ASSEMBLY_DIR/final.contigs.fa -o $ASSEMBLY_DIR/contigs.fa --min-len $MIN_CONTIG_SIZE --simplify-names --report $ASSEMBLY_DIR/name_conversions.txt

# gzip $ASSEMBLY_DIR/final.contigs.fa

bowtie2-build $ASSEMBLY_DIR/contigs.fa $ASSEMBLY_DIR/contigs
bowtie2 --threads $NUM_THREADS -x $ASSEMBLY_DIR/contigs -1 $R1s -2 $R2s -S $ASSEMBLY_DIR/$ASB.sam

cd $ASSEMBLY_DIR

ls /users/work/cat3/projects/mime/sunbeam_output/qc/decontam/*_wgs_1.fastq.gz | grep "wgs" | grep -e 1045m -e 1004m -e 90m -e 470m -e MC | cut -f 10 -d "/" | cut -f 1,2,3,4,5,6 -d "_" > samples.txt

for sample in `awk '{print $1}' samples.txt`
do
    if [ "$sample" == "sample" ]; then continue; fi
    
    # you need to make sure you "ls 01_QC/*QUALITY_PASSED_R1*" returns R1 files for all your samples in samples.txt
    R1s=$(ls /users/work/cat3/projects/mime/sunbeam_output/qc/decontam/$sample*1.fastq.gz | python -c 'import sys; print( ",".join([x.strip() for x in sys.stdin.readlines()]))' )
    R2s=$(ls /users/work/cat3/projects/mime/sunbeam_output/qc/decontam/$sample*2.fastq.gz | python -c 'import sys; print( ",".join([x.strip() for x in sys.stdin.readlines()]))' )

    bowtie2 --threads $NUM_THREADS -x $ASSEMBLY_DIR/contigs -1 $R1s -2 $R2s --no-unal -S $ASSEMBLY_DIR/$sample.sam
    samtools view -F 4 -bS $ASSEMBLY_DIR/$sample.sam > $ASSEMBLY_DIR/$sample-RAW.bam
    anvi-init-bam $ASSEMBLY_DIR/$sample-RAW.bam -o $ASSEMBLY_DIR/$sample.bam
    rm $ASSEMBLY_DIR/$sample.sam $ASSEMBLY_DIR/$sample-RAW.bam
done

# mapping is done, and we no longer need bowtie2-build files
rm $ASSEMBLY_DIR/*.bt2

anvi-gen-contigs-database -f $ASSEMBLY_DIR/contigs.fa -o $ASSEMBLY_DIR/contigs.db -n $ASB

anvi-display-contigs-stats $ASSEMBLY_DIR/contigs.db --report-as-text --output-file contigs-stats.txt

anvi-run-hmms -c $ASSEMBLY_DIR/contigs.db -T $NUM_THREADS

anvi-run-ncbi-cogs -c $ASSEMBLY_DIR/contigs.db --num-threads $NUM_THREADS

for sample in `cat samples.txt`; do 
  anvi-profile -i $ASSEMBLY_DIR/$sample.bam -c $ASSEMBLY_DIR/contigs.db -T $NUM_THREADS -o $PROFILE_DIR/$sample --sample-name "Sample_"$sample --profile-SCVs --cluster-contigs
done

anvi-merge $PROFILE_DIR/*/PROFILE.db -o $MERGE_DIR -c $ASSEMBLY_DIR/contigs.db --enforce-hierarchical-clustering -S $ASB -W

# ----- TAXONOMY ANNOTATIONS KAIJU -----

anvi-get-sequences-for-gene-calls -c $ASSEMBLY_DIR/contigs.db -o $ASSEMBLY_DIR/gene_calls.fa

kaiju -t /users/work/cat3/db/kaiju/nodes.dmp \
      -f /users/work/cat3/db/kaiju/kaiju_db_nr_euk.fmi \
      -i $ASSEMBLY_DIR/gene_calls.fa \
      -o $ASSEMBLY_DIR/gene_calls.nr_euk.out \
      -z $NUM_THREADS \
      -v

kaiju-addTaxonNames -t /users/work/cat3/db/kaiju/nodes.dmp \
              -n /users/work/cat3/db/kaiju/names.dmp \
              -i $ASSEMBLY_DIR/gene_calls.nr_euk.out \
              -o $ASSEMBLY_DIR/gene_calls.nr_euk.names \
              -r superkingdom,phylum,order,class,family,genus,species

anvi-import-taxonomy-for-genes -i $ASSEMBLY_DIR/gene_calls.nr_euk.names -c $ASSEMBLY_DIR/contigs.db -p kaiju --just-do-it

# ----- BINNING -----

# concoct
anvi-cluster-contigs -p $MERGE_DIR/PROFILE.db -c $ASSEMBLY_DIR/contigs.db -C "concoct" --driver concoct -T $NUM_THREADS --just-do-it
anvi-estimate-genome-completeness -c $ASSEMBLY_DIR/contigs.db -p $MERGE_DIR/PROFILE.db -C "concoct" -o $ASSEMBLY_DIR/estimate-genome-completeness

# maxbin2
anvi-cluster-contigs -p $MERGE_DIR/PROFILE.db -c $ASSEMBLY_DIR/contigs.db -C "maxbin2" --driver maxbin2 -T $NUM_THREADS --just-do-it
anvi-estimate-genome-completeness -c $ASSEMBLY_DIR/contigs.db -p $MERGE_DIR/PROFILE.db -C "maxbin2" -o $ASSEMBLY_DIR/estimate-genome-completeness

# metabat2
anvi-cluster-contigs -p $MERGE_DIR/PROFILE.db -c $ASSEMBLY_DIR/contigs.db -C "metabat2" --driver metabat2 -T $NUM_THREADS --just-do-it
anvi-estimate-genome-completeness -c $ASSEMBLY_DIR/contigs.db -p $MERGE_DIR/PROFILE.db -C "metabat2" -o $ASSEMBLY_DIR/estimate-genome-completeness

# ----- SUMMARY -----

anvi-summarize -p $MERGE_DIR/PROFILE.db \
    -c $ASSEMBLY_DIR/contigs.db \
    -C "concoct" \
    -o $SUMMARY_DIR/concoct

anvi-summarize -p $MERGE_DIR/PROFILE.db \
    -c $ASSEMBLY_DIR/contigs.db \
    -C "maxbin2" \
    -o $SUMMARY_DIR/maxbin2

anvi-summarize -p $MERGE_DIR/PROFILE.db \
    -c $ASSEMBLY_DIR/contigs.db \
    -C "metabat2" \
    -o $SUMMARY_DIR/metabat2

# ----- PROKKA ANNOTATIONS -----

prokka --prefix PROKKA \
       --outdir $ASSEMBLY_DIR/prokka \
       --cpus $NUM_THREADS \
       --metagenome $ASSEMBLY_DIR/contigs.fa

python /users/work/cat3/scripts/pygff_parser.py $ASSEMBLY_DIR/prokka/PROKKA.gff \
                         --gene-calls $ASSEMBLY_DIR/prokka/gene_calls.txt \
                         --annotation $ASSEMBLY_DIR/prokka/gene_annot.txt

anvi-gen-contigs-database -f $ASSEMBLY_DIR/contigs.fa \
  -o $ASSEMBLY_DIR/contigs-prokka.db \
  --external-gene-calls $ASSEMBLY_DIR/prokka/gene_calls.txt \
  --project-name 'WGS bottom seawater' \
  --ignore-internal-stop-codons

anvi-import-functions -c $ASSEMBLY_DIR/contigs.db \
  -i $ASSEMBLY_DIR/prokka/gene_annot.txt

# ----- Pfam FUNCTIONS ------

anvi-run-pfams -c $ASSEMBLY_DIR/contigs.db -T $NUM_THREADS

anvi-export-functions -c $ASSEMBLY_DIR/contigs.db --list-annotation-sources

anvi-interactive -p $MERGE_DIR/PROFILE.db -c $ASSEMBLY_DIR/contigs.db --server-only -P 8080

# ---- FILTER MAGS FROM BINS -----

BACTDIR=/users/work/cat3/projects/mime/results/de_novo_wf/bottom-bact-mags
mkdir -p $BACTDIR

for id in `cat $SUMMARY_DIR/metabat2/bins_summary.txt | sed -e 's/ /_/g' | awk '{if (($7 >= 50) && ($8 <= 10)) { print } }' | grep METABAT | cut -f 1` ; do
  domain=$(head "$SUMMARY_DIR/metabat2/bin_by_bin/$id/$id-scg_domain.txt")
  if [[ "$domain" = bacteria ]]; then
    echo $id $domain
    cp $SUMMARY_DIR/metabat2/bin_by_bin/$id/$id-contigs.fa $BACTDIR
  fi
done

ARCHDIR=/users/work/cat3/projects/mime/results/de_novo_wf/bottom-arch-mags
mkdir -p $ARCHDIR

archaea_count=0
bacteria_count=0

for id in `cat $SUMMARY_DIR/metabat2/bins_summary.txt | sed -e 's/ /_/g' | awk '{if (($7 >= 50) && ($8 <= 10)) { print } }' | grep METABAT | cut -f 1` ; do
  domain=$(head "$SUMMARY_DIR/metabat2/bin_by_bin/$id/$id-scg_domain.txt")
  if [[ "$domain" = archaea ]]; then
    echo $id $domain
    archaea_count=$((archaea_count+1))
    cp $SUMMARY_DIR/metabat2/bin_by_bin/$id/$id-contigs.fa $ARCHDIR
  elif [[ "$domain" = bacteria ]]; then
    bacteria_count=$((bacteria_count+1))
  fi
done

echo $archaea_count
echo $bacteria_count

# ready to use GTDB-Tk :)
