# simkaviz_pca.r

# ------ PRELIMINARY SETTINGS ------

# import colours to remain consistent between all plots
source("/Users/Clara/Projects/colors/colors.R")
source("/Users/Clara/Projects/colors/colors2.R")

# set working directory
setwd("~/mime/results/simkamin")

library(vegan)
library(ggvegan)
library(ggplot2)
library(ggpubr)

# set the directory to save img
img.path = "~/graphics/plots/"
library(stringr)
today <- str_replace_all(Sys.Date(), "-", "") 

# ---- IMPORT DATA -----

# Distance matrix
distanceMatrixFilename = "mat_abundance_braycurtis.csv.gz"

# Metadata
metadata_table = as.matrix(read.table(file="metadata_simka.csv", sep=";", header=TRUE, row.names=1))

# ---- CHOOSE PCA AXIS -----

pca_axis1 = 1
pca_axis2 = 2

# ----- READ MATRIX -----

distanceMatrix = as.matrix(read.table(file=distanceMatrixFilename, sep=";", header=TRUE, row.names=1))
distanceMatrix[lower.tri(distanceMatrix)] <- t(distanceMatrix)[lower.tri(distanceMatrix)] #symmetrize matrix

distance_name = basename(distanceMatrixFilename)
distance_name = unlist(strsplit(distance_name, "[.]"))[1]
distance_name = gsub("mat_", "", distance_name)

pcoa <- cmdscale(as.dist(distanceMatrix), k = 2, eig = TRUE)

df <- as.data.frame(cbind(pcoa$points[,1], pcoa$points[,2]))
colnames(df) <- c("x", "y")
df$date <- as.Date(metadata_table[,c("Date")],"%Y%m%d")
df$stn <- metadata_table[,c("Station")]
df$reg <- metadata_table[,c("Region")]
df$run <- metadata_table[,c("Run")]
df$lib.prep <- metadata_table[,c("Lib")]
df$depth <- metadata_table[,c("Depth")]
df$cat <- "seafloor"
df[df$depth == "0m",]$cat <- "surface"

getSeason <- function(input.date){
  numeric.date <- 100*month(input.date)+day(input.date)
  ## input Seasons upper limits in the form MMDD in the "break =" option:
  cuts <- base::cut(numeric.date, breaks = c(0,319,0620,0921,1220,1231)) 
  # rename the resulting groups (could've been done within cut(...levels=) if "Winter" wasn't double
  levels(cuts) <- c("Winter","Spring","Summer","Fall","Winter")
  return(cuts)
}

df <- df %>% as_tibble() %>%
  mutate(season = getSeason(date))
        

ggplot(df, aes(x=df$x, y=df$y)) + 
  geom_point(data = df, aes(col = season, shape=lib.prep)) +
  #geom_text(aes(label=df$stn, col=df$stn), hjust=0,vjust=0) +
  theme(aspect.ratio=1) 
  #theme_pubr()

ggplot(df, aes(x=df$x, y=df$y)) + 
  geom_point(aes(col=ceiling_date(df$date, "season"), shape=df$lib)) +
  #geom_text(aes(label=df$stn, col=df$stn), hjust=0,vjust=0) +
  theme(aspect.ratio=1) 
#theme_pubr()

adonis(distanceMatrix ~ df$stn, perm=200)
adonis(distanceMatrix ~ metadata_table[,c("Depth")], perm=200)
adonis(distanceMatrix ~ metadata_table[,c("Station")], perm=200)
adonis(distanceMatrix ~ metadata_table[,c("Run")], perm=200)
adonis(distanceMatrix ~ metadata_table[,c("Lib")], perm=200)
adonis(distanceMatrix ~ metadata_table[,c("Depth")], perm=200)

# ----- WITHOUT MOCK -----

distanceMatrix2 <- broom::tidy(as.dist(distanceMatrix)) %>% 
  as_tibble() %>%
  filter(item1 != "ID32" & item2 != "ID37") %>%
  tidyr::spread(item1, value = distance, fill = 0)

mock <- c("ID32", "ID37")
samples <- setdiff(colnames(distanceMatrix), mock)
metadata_table3 <- metadata_table[samples,] %>% as.data.frame() %>% droplevels()

distanceMatrix3 <- distanceMatrix[samples, samples]
# distanceMatrix3 <- as.matrix(distanceMatrix2[,!colnames(distanceMatrix2) %in% c("ID32", "ID37")])
# distanceMatrix3 <- apply(distanceMatrix3, 1, as.numeric) 

adonis(distanceMatrix3 ~ metadata_table3[,c("Region")], perm=200)
adonis(distanceMatrix3 ~ metadata_table3[,c("Depth")], perm=200)
adonis(distanceMatrix3 ~ metadata_table3[,c("Station")], perm=200)
adonis(distanceMatrix3 ~ metadata_table3[,c("Run")], perm=200)
adonis(distanceMatrix3 ~ metadata_table3[,c("Lib")], perm=200)
adonis(distanceMatrix3 ~ metadata_table3[,c("Depth")], perm=200)

pcoa <- cmdscale(as.dist(as.data.frame(distanceMatrix3)), k = 2, eig = TRUE)

df <- as.data.frame(cbind(pcoa$points[,1], pcoa$points[,2]))
colnames(df) <- c("x", "y")
df$stn <- metadata_table3[,c("Station")]
df$reg <- metadata_table3[,c("Region")]
df$run <- metadata_table3[,c("Run")]
df$lib.prep <- metadata_table3[,c("Lib")]
df$depth <- metadata_table3[,c("Depth")]
df$cat <- "seafloor"
df[df$depth == "0m",]$cat <- "surface"

lab.plot <- ggplot(df, aes(x=df$x, y=df$y)) + 
  geom_point(size=5, aes(col=df$run, shape=df$lib)) +
  #geom_text(aes(label=df$stn), col="black", hjust=.5,vjust=.5) +
  theme_pubr() + theme(aspect.ratio=1, legend.position="right") + 
  scale_color_manual(values=c(Sunflower, Mint)) +
  labs(title = "SimkaMin PCoA", subtitle = "Bray curtis distances", 
       col = "HiSeq run", shape="Library prep.",
       x="PC1", y="PC2")
lab.plot


pdf(file = paste(img.path, today, "_simka_min_lab.pdf", sep="", collapse=NULL),  
    width=8, height=6)
lab.plot
dev.off()

env.plot <- ggplot(df, aes(x=df$x, y=df$y)) + 
  geom_point(size=5, aes(col=df$stn, shape=df$reg)) +
  geom_text(aes(label=df$stn), col="black", hjust=.5,vjust=.5) +
  theme_pubr() + theme(aspect.ratio=1, legend.position="right") +  
  scale_color_manual(values=c(Grapefruit, Orange, Jeans, Aqua)) +
  #scale_shape_manual(values=c(78, 83)) +
  labs(title = "SimkaMin PCoA", subtitle = "Bray curtis distances",
       col = "Station ID", shape="Region", 
       x="PC1", y="PC2")
env.plot

# UPDATE
env.plot <- ggplot(df, aes(x=df$x, y=df$y)) + 
  geom_point(size=5, aes(col=df$stn, shape=df$cat)) +
  geom_text(aes(label=df$stn), col="white", hjust=.5,vjust=.5, size=1.5) +
  theme_pubr() + theme(aspect.ratio=1, legend.position="right") +  
  scale_color_manual(values=c(Grapefruit, Orange, Jeans, Aqua)) +
  #scale_shape_manual(values=c(78, 83)) +
  labs(title = "SimkaMin PCoA", subtitle = "Bray curtis distances",
       col = "Station ID", shape="Environment", 
       x="PC1", y="PC2")
env.plot

# FOCUS ON WGS
df <- df[df$lib == "wgs",]

lab.plot <- ggplot(df, aes(x=df$x, y=df$y)) + 
  geom_point(size=5, stroke = 1, shape = 21, aes(fill=run)) +
  geom_text(aes(label=df$stn), col="white", hjust=.5,vjust=.5, size=1.5) +
  scale_fill_manual(values=c(Sunflower, Mint)) +
  #scale_shape_manual(values=c(78, 83)) +
  theme_void() +
  labs(title = "SimkaMin PCoA", subtitle = "Bray curtis distances",
       fill = "HiSeq Run", 
       x="PC1", y="PC2") +
  theme(aspect.ratio=1,
        panel.border = element_rect(colour = "black", fill=NA, size=1.5),
        axis.ticks = element_line(size=.8, color="black"),
        axis.ticks.length=unit(0.2,"cm"),
        legend.box = "vertical",
        legend.position = c(0.85, 0.15),
        text=element_text(family="Frutiger"))  + font("xylab",size=20, face = "bold") +
  font("xy", size=12, face = "bold") +
  font("xy.text", size = 12) +
  font("legend.title",size = 10, face = "bold") +
  font("legend.text",size = 10)
lab.plot

env.plot <- ggplot(df, aes(x=df$x, y=df$y)) + 
  geom_point(size=5, stroke = 1, aes(fill=cat, shape=reg)) +
  geom_text(aes(label=df$stn), col="white", hjust=.5,vjust=.5, size=1.5) +
  scale_shape_manual(values=c(21, 23)) +
  scale_fill_manual(values=c(Jeans, Grass)) +
  #scale_shape_manual(values=c(78, 83)) +
  theme_void() +
  labs(title = "SimkaMin PCoA", subtitle = "Bray curtis distances",
       fill = "Environment", shape="Region", 
       x="PC1", y="PC2") +
  theme(aspect.ratio=1, legend.position = c(0.7, 0.15),
        panel.border = element_rect(colour = "black", fill=NA, size=1.5),
        axis.ticks = element_line(size=.8, color="black"),
        axis.ticks.length=unit(0.2,"cm"),
        legend.box = "horizontal",
        text=element_text(family="Frutiger"))  + font("xylab",size=20, face = "bold") +
  font("xy", size=12, face = "bold") +
  font("xy.text", size = 12) +
  font("legend.title",size = 10, face = "bold") +
  font("legend.text",size = 10)
env.plot

plot_grid(lab.plot, env.plot, labels = "AUTO")

pdf(file = paste(img.path, today, "_simka_min_env.pdf", sep="", collapse=NULL),  
    width=8, height=6)
env.plot
dev.off()

