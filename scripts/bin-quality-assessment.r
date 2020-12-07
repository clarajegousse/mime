# bin-quality-assessment.r

# ----- PRELIMINARY SETTINGS -----

library(devtools)
# import colours to remain consistent between all plots
source_url("https://raw.githubusercontent.com/clarajegousse/colors/master/colors.R")

# set the directory to save img
img.path = "/Users/Clara/Projects/diary/graphics/plots/"

library(stringr)
today <- str_replace_all(Sys.Date(), "-", "") # concatenate the date following the format "yyyymmdd"
# the date is used to save output files such as plots without overwritting previous work (on a daily basis)

# ----- LIBRARIES -----

library(ggpubr)
library(dplyr)
library(ggrepel)
library(ggplot2)
library(cowplot)
library(extrafont)

# ----- TRANSFER DATA FROM GARPUR TO LOCAL -----

# ------ SURFACE ---------

# ------ IMPORT DATA -----

file.path <- "https://raw.githubusercontent.com/clarajegousse/mime/main/results/coassembly/coassembly_wgs_surface/summary/metabat2/bins_summary.txt"
df <- read.table(file.path, 
                 header = TRUE, 
                 sep = "\t", 
                 dec = ".")
head(df)
df$taxon <- as.character(df$taxon)
sapply(df, mode)

# quick overview
ggscatter(df, x = "percent_completion", y = "percent_redundancy")

# ----- BIN TYPES -----

df$bin.type <- "bins"
df[df$percent_completion >= 50 & df$percent_redundancy <= 10,]$bin.type <- "mags"

table(df$bin.type)

# ----- BIN QUALITY CATEGORIES ----

df$quality <- "bad"

df[df$percent_completion >= 50 & df$percent_redundancy <= 10,]$quality <- "good"

# df[df$percent_completion >= 50 & df$percent_redundancy <= 4,]$quality <- "partial"
# df[df$percent_completion >= 70 & df$percent_redundancy <= 10,]$quality <- "medium"
# df[df$percent_completion >= 90 & df$percent_redundancy <= 5,]$quality <- "near complete"

df$quality <- as.factor(df$quality)
df$quality <- ordered(df$quality)

# df$quality <- ordered(df$quality, levels = c("bad", "partial", "medium", "near complete"))
df$quality <- ordered(df$quality, levels = c("bad", "good"))


table(df$quality)

lvls = levels(as.factor(df$quality))
labels = paste(lvls," (",table(df$quality),")",sep="")
# col.palette <- c(HIred, HIorange, HIyellow, HIgreen)
# col.palette <- c(Grapefruit, Orange, Sunflower, Grass)
col.palette <- c(Grapefruit, Mint)

df.surface <- df

# ----- HISTOGRAM -----

pmain <- ggplot(df, aes(x = percent_completion, 
                        y = percent_redundancy, 
                        fill = quality, 
                        size = total_length)) +
  geom_point(aes(fill = quality), shape = 21, stroke = .5) + theme_pubr() + theme(aspect.ratio=1) +
  scale_fill_manual(values = col.palette, labels = labels) +
  ylab(c("Redundancy (%)")) +
  xlab(c("Completion (%)")) +
  labs(fill= "Quality", size = "Total length (pb)") +
  #rremove("legend.title") +
  theme(panel.border = element_rect(colour = "black", fill=NA, size=.8),
        axis.line = element_line(size=0,color="red"),
        axis.ticks = element_line(size=.5,color="black"),
        axis.ticks.length=unit(0.2,"cm"),
        #legend.position = c(0.2, 0.7),
        legend.position = "right",
        text=element_text(family="Frutiger"))  + font("xylab",size = 20, face = "bold") +
  font("xy", size=12, face = "bold") +
  font("xy.text", size = 12) +
  font("legend.title",size = 12, face = "bold") +
  font("legend.text",size = 10) + 
  guides(fill = guide_legend(override.aes = list(size = 5)))+ 
  annotate("rect", xmin = 49, xmax = 101, ymin = -1, ymax = 10, 
           linetype = 1, size = .5, color = "black", alpha = 0) + ylim(0, 101)


xhist <- axis_canvas(pmain, axis = "x")+
  geom_histogram(data = df, aes(x = percent_completion, fill = quality),
                 size = .5, color = "black") +
  scale_fill_manual(values = col.palette)

yhist <- axis_canvas(pmain, axis = "y", coord_flip = TRUE)+
  geom_histogram(data = df, aes(x = percent_redundancy, fill = quality),
                size = .5, color = "black") +
  scale_fill_manual(values = col.palette) +
  coord_flip()

p1 <- insert_xaxis_grob(pmain, xhist, grid::unit(.2, "null"), position = "top")
p2 <- insert_yaxis_grob(p1, yhist, grid::unit(.2, "null"), position = "right")

final <- ggdraw(p2)
final 

plot1.surface <- final
plot1.surface
# save as svg heigh=600

# file.format <- "png"
# file.name <- paste(img.path, today, "_coassembly_wgs_surface_assessment_of_bins_quality_metabat2.", file.format, sep="", collapse=NULL)
# ggsave(final, file=file.name, dpi = 300, units = "cm", width = 22)

# ----- FOCUS ON MOCK -----

df.mags <- df[df$percent_completion >= 50 & df$percent_redundancy <= 10,]

table(df.mags$quality)
lvls = levels(as.factor(df.mags$quality))
labels = paste(lvls," (",table(df.mags$quality),")",sep="")

mock.genus <- c("Sulfitobacter", 
             "Alteromonas", 
             "Bacillus",
             "Geoacillus",
             "Colwellia",
             "Dietzia",
             "Escherichia",
             "Halomonas",
             "Marinobacter",
             "Photobacterium",
             "Pseudoalteromonas",
             "Pyrococcus",
             "Reinekea",
             "Rhodococcus",
             "Thermococcus",
             "Thermus",
             "Vibrio")

df.mags$community <- NA
df.mags[df.mags$taxon %in% mock.genus,]$community <- "mock"

`%notin%` <- Negate(`%in%`)
df.mags[df.mags$taxon %notin% mock.genus,]$community <- "natural"

# corrections
# METABAT__176 is clearly not in the mock community according to differential coverage
df.mags[df.mags$bins == "METABAT__176",]$community <- "natural"

df.mags[df.mags$bins == "METABAT__221",]$community <- "mock" # Escherichia coli
df.mags[df.mags$bins == "METABAT__221",]$taxon <- "Escherichia"

df.mags[df.mags$bins == "METABAT__203",]$community <- "mock" # Geobacillus
df.mags[df.mags$bins == "METABAT__140",]$community <- "mock" # Pseudomonas
table(df.mags$community)
lvls2 = levels(as.factor(df.mags$community))
labels2 = paste(lvls2," (",table(df.mags$community),")",sep="")


pmain <- ggplot(df.mags, aes(x = percent_completion, 
                        y = percent_redundancy, 
                        label = taxon,
                        #fill = quality, 
                        size = total_length)) +
  geom_point(aes(fill = community), shape = 21, stroke = 1) + theme_pubr() + theme(aspect.ratio=1) +
  scale_fill_manual(values = col.palette, labels = labels2) +
  #scale_shape_manual(values=c(21, 23), labels = labels2) +
  ylab(c("Redundancy (%)")) +
  xlab(c("Completion (%)")) +
  labs(fill= "Community", size = "Total length (pb)") +
  #rremove("legend.title") +
  theme(panel.border = element_rect(colour = "black", fill=NA, size=1.5),
        axis.ticks = element_line(size=.8, color="black"),
        axis.ticks.length=unit(0.2,"cm"),
        legend.box = "vertical",
        legend.position = "right", #c(0.2, 0.7),
        text=element_text(family="Frutiger"))  + font("xylab",size=20, face = "bold") +
  font("xy", size=12, face = "bold") +
  font("xy.text", size = 12) +
  font("legend.title",size = 10, face = "bold") +
  font("legend.text",size = 10)

plot2.surface <- pmain + xlim(48, 100) + ylim(0, 10) +
  geom_text_repel(
    data = subset(df.mags, community == "mock"),
    size = 4,
    max.iter = 2000, force =100, point.padding = 1 ) +
  guides(fill = guide_legend(override.aes = list(size = 4, fill = col.palette, shape = 21, stroke = 1)),
         size = guide_legend(override.aes = list(shape = 21)))

plot2.surface

# ------ SEA FLOOR ---------

# ------ IMPORT DATA -----

file.path <- "~/Projects/mime/results/coassembly/coassembly_wgs_bottom/bins_summary.txt"
df <- read.table(file.path, 
                 header = TRUE, 
                 sep = "\t", 
                 dec = ".")
head(df)
df$taxon <- as.character(df$taxon)
sapply(df, mode)

# quick overview
ggscatter(df, x = "percent_completion", y = "percent_redundancy")

# ----- BIN TYPES -----

df$bin.type <- "bins"
df[df$percent_completion >= 50 & df$percent_redundancy <= 10,]$bin.type <- "mags"

table(df$bin.type)

# ----- BIN QUALITY CATEGORIES ----

df$quality <- "bad"

df[df$percent_completion >= 50 & df$percent_redundancy <= 10,]$quality <- "good"

#df[df$percent_completion >= 50 & df$percent_redundancy <= 4,]$quality <- "partial"
#df[df$percent_completion >= 70 & df$percent_redundancy <= 10,]$quality <- "medium"
#df[df$percent_completion >= 90 & df$percent_redundancy <= 5,]$quality <- "near complete"

df$quality <- as.factor(df$quality)
df$quality <- ordered(df$quality)
# df$quality <- ordered(df$quality, levels = c("bad", "partial", "medium", "near complete"))

df$quality <- ordered(df$quality, levels = c("bad", "good"))


table(df$quality)

lvls = levels(as.factor(df$quality))
labels = paste(lvls," (",table(df$quality),")",sep="")
# col.palette <- c(HIred, HIorange, HIyellow, HIgreen)
# col.palette <- c(Grapefruit, Orange, Sunflower, Grass)
col.palette <- c(Grapefruit, Mint)

df.seafloor <- df

# ----- HISTOGRAM -----

pmain <- ggplot(df, aes(x = percent_completion, 
                        y = percent_redundancy, 
                        fill = quality, 
                        size = total_length)) +
  geom_point(aes(fill = quality), shape = 21, stroke = .5) + theme_pubr() + theme(aspect.ratio=1) +
  scale_fill_manual(values = col.palette, labels = labels) +
  ylab(c("Redundancy (%)")) +
  xlab(c("Completion (%)")) +
  labs(fill= "Quality", size = "Total length (pb)") +
  #rremove("legend.title") +
  theme(panel.border = element_rect(colour = "black", fill=NA, size=.8),
        axis.line = element_line(size=0,color="red"),
        axis.ticks = element_line(size=.5,color="black"),
        axis.ticks.length=unit(0.2,"cm"),
        legend.position = c(0.2, 0.7),
        text=element_text(family="Frutiger"))  + font("xylab",size = 20, face = "bold") +
  font("xy", size=12, face = "bold") +
  font("xy.text", size = 12) +
  font("legend.title",size = 12, face = "bold") +
  font("legend.text",size = 10) + 
  guides(fill = guide_legend(override.aes = list(size = 5)))+ 
  annotate("rect", xmin = 49, xmax = 101, ymin = -1, ymax = 10, 
           linetype = 1, size = .5, color = "black", alpha = 0)

pmain + ylim(0, 100)

xhist <- axis_canvas(pmain, axis = "x")+
  geom_histogram(data = df, aes(x = percent_completion, fill = quality),
                 size = .5, color = "black") +
  scale_fill_manual(values = col.palette)

yhist <- axis_canvas(pmain, axis = "y", coord_flip = TRUE)+
  geom_histogram(data = df, aes(x = percent_redundancy, fill = quality),
                 size = .5, color = "black") +
  scale_fill_manual(values = col.palette) +
  coord_flip()

p1 <- insert_xaxis_grob(pmain, xhist, grid::unit(.2, "null"), position = "top")
p2 <- insert_yaxis_grob(p1, yhist, grid::unit(.2, "null"), position = "right")

final <- ggdraw(p2)
final 

plot1.bottom <- final
plot1.bottom

# file.format <- "png"
# file.name <- paste(img.path, today, "_coassembly_wgs_surface_assessment_of_bins_quality_metabat2.", file.format, sep="", collapse=NULL)
# ggsave(final, file=file.name, dpi = 300, units = "cm", width = 22)

# ----- FOCUS ON MOCK -----

df.mags <- df[df$percent_completion >= 50 & df$percent_redundancy <= 10,]

table(df.mags$quality)
lvls = levels(as.factor(df.mags$quality))
labels = paste(lvls," (",table(df.mags$quality),")",sep="")

mock.genus <- c("Sulfitobacter", 
                "Alteromonas", 
                "Bacillus",
                "Geobacillus",
                "Colwellia",
                "Dietzia",
                "Escherichia",
                "Halomonas",
                "Marinobacter",
                "Photobacterium",
                "Pseudoalteromonas",
                "Pyrococcus",
                "Reinekea",
                "Rhodococcus",
                "Thermococcus",
                "Thermus",
                "Vibrio")

df.mags$community <- NA
df.mags[df.mags$taxon %in% mock.genus,]$community <- "mock"

`%notin%` <- Negate(`%in%`)
df.mags[df.mags$taxon %notin% mock.genus,]$community <- "natural"


# corrections
# METABAT__176 is clearly not in the mock community according to differential coverage
df.mags[df.mags$bins == "METABAT__138",]$community <- "mock" # Eschierichia coli
df.mags[df.mags$bins == "METABAT__138",]$taxon <- "Escherichia"

table(df.mags$community)
lvls2 = levels(as.factor(df.mags$community))
labels2 = paste(lvls2," (",table(df.mags$community),")",sep="")


pmain <- ggplot(df.mags, aes(x = percent_completion, 
                             y = percent_redundancy, 
                             label = taxon,
                             #fill = quality, 
                             size = total_length)) +
  geom_point(aes(fill = community), shape = 21, stroke = 1) + theme_pubr() + theme(aspect.ratio=1) +
  scale_fill_manual(values = col.palette, labels = labels2) +
  #scale_shape_manual(values=c(21, 23), labels = labels2) +
  ylab(c("Redundancy (%)")) +
  xlab(c("Completion (%)")) +
  labs(fill= "Community", size = "Total length (pb)") +
  #rremove("legend.title") +
  theme(panel.border = element_rect(colour = "black", fill=NA, size=1.5),
        axis.ticks = element_line(size=.8, color="black"),
        axis.ticks.length=unit(0.2,"cm"),
        legend.box = "vertical",
        legend.position = "right", #c(0.2, 0.7),
        text=element_text(family="Frutiger"))  + font("xylab",size=20, face = "bold") +
  font("xy", size=12, face = "bold") +
  font("xy.text", size = 12) +
  font("legend.title",size = 10, face = "bold") +
  font("legend.text",size = 10)

plot2.bottom <- pmain + xlim(48, 100) + ylim(0, 10) +
  geom_text_repel(
    data = subset(df.mags, community == "mock"),
    size = 4,
    max.iter = 2000, force =100, point.padding = 1 ) +
  guides(fill = guide_legend(override.aes = list(size = 4, fill = col.palette, shape = 21, stroke = 1)),
         size = guide_legend(override.aes = list(shape = 21)))

plot2.bottom

# ----- FINAL PLOT -----

plot_grid(plot1.surface, plot1.bottom,
  labels = "AUTO")


plot_grid(plot2.surface, plot2.bottom,
          labels = "AUTO")

# ----- CHI SQUARE TEST -----

# prep of the global dataset
df.surface$cat <- "surface"

df.surface$community <- NA
df.surface[df.surface$taxon %in% mock.genus,]$community <- "mock"

`%notin%` <- Negate(`%in%`)
df.surface[df.surface$taxon %notin% mock.genus,]$community <- "natural"

# corrections
# METABAT__176 is clearly not in the mock community according to differential coverage
df.surface[df.surface$bins == "METABAT__176",]$community <- "natural"

df.surface[df.surface$bins == "METABAT__221",]$community <- "mock" # Escherichia coli
df.surface[df.surface$bins == "METABAT__221",]$taxon <- "Escherichia"

df.surface[df.surface$bins == "METABAT__203",]$community <- "mock" # Geobacillus
df.surface[df.surface$bins == "METABAT__140",]$community <- "mock" # Pseudomonas

df.seafloor$cat <- "seafloor"

df.seafloor$community <- NA
df.seafloor[df.seafloor$taxon %in% mock.genus,]$community <- "mock"

`%notin%` <- Negate(`%in%`)
df.seafloor[df.seafloor$taxon %notin% mock.genus,]$community <- "natural"

# corrections
# METABAT__176 is clearly not in the mock community according to differential coverage
df.seafloor[df.seafloor$bins == "METABAT__138",]$community <- "mock" # Eschierichia coli
df.seafloor[df.seafloor$bins == "METABAT__138",]$taxon <- "Escherichia"

df <- rbind.data.frame(df.surface, df.seafloor)

df$community <- NA
df[df$taxon %in% mock.genus,]$community <- "mock"

`%notin%` <- Negate(`%in%`)
df[df$taxon %notin% mock.genus,]$community <- "natural"

# chi square test

# H0: The two variables are independent // the quality of the bins and the env. are independant
# H1: The two variables relate to each other // the quality of the bins depends on the environment

table(df$quality, df$cat)
chisq.test(df$quality, df$cat, correct=TRUE)

# We have a chi-squared value of 0.27784 and 
# since we get a p-Value above than the significance level of 0.05, 
# we cannot reject the null hypothesis and conclude that the two variables are in fact dependant

table(df[df$quality=="good",]$community, df[df$quality=="good",]$cat)
chisq.test(df[df$quality=="good",]$community, df[df$quality=="good",]$cat, correct=TRUE)
