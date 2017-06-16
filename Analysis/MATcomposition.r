# Paul Selmants
# December 21, 2013, modified March 31, 2015
## MAT bacteria ordinations and permanova ##

#load required packages
library(vegan)
library(tidyr)
library(ggplot2)
library(dplyr)

#read distance matrices into R (Bray, weighted UniFrac and unweighted UniFrac)
bray.dm <- read.delim('bray_curtis_dm.txt', row.names = 1)
wUF.dm <- read.delim('weighted_unifrac_dm.txt', row.names = 1)

#perform global NMDS ordinations using monoMDS with bray and wUNIFRAC distance matrices
bray.ord <- metaMDS(bray.dm)
wUF.ord <- metaMDS(wUF.dm)
 
#read environmental matrix into R 
site.env <- read.csv('site_mat_sorted.csv', row.names = 1)
#test effect of MAT on community composition using permanova 
bray.mod <- adonis(bray.dm ~ MAT + TBCF + pH, site.env, permutations = 999)
wUF.mod <- adonis(wUF.dm ~ MAT + TBCF + pH, site.env, permutations = 999)

#create dataframe including bray NMDS axis scores & site data 
bray.xy <- data.frame(MDS1 = bray.ord$points[,1], MDS2 = bray.ord$points[,2])
bray.xy$sampleID = rownames(bray.xy)
rownames(bray.xy) = NULL
bray.MAT <- bray.xy %>%
  separate(sampleID, into = c('plot', 'rep'), sep = '\\.') %>%
  mutate(MAT = c(rep(13.8, 8), rep(13.0, 6), rep(16.7, 8), rep(16.1, 8), 
    rep(18.2, 8), rep(17.3, 8),rep(16.1, 8), rep(15.5, 8), rep(15.1, 7))) %>%
  arrange(MAT, plot, rep) %>%
  select(plot, rep, MAT, MDS1, MDS2)

#create dataframe including weighted Unifrac NMDS axis scores & site data
wUF.xy <- data.frame(MDS1 = wUF.ord$points[,1], MDS2 = wUF.ord$points[,2])
wUF.xy$sampleID = rownames(wUF.xy)
rownames(wUF.xy) = NULL
wUF.MAT <- wUF.xy %>%
  separate(sampleID, into = c('plot', 'rep'), sep = '\\.') %>%
  mutate(MAT = c(rep(13.8, 8), rep(13.0, 6), rep(16.7, 8), rep(16.1, 8), 
    rep(18.2, 8), rep(17.3, 8),rep(16.1, 8), rep(15.5, 8), rep(15.1, 7))) %>%
  arrange(MAT, plot, rep) %>%
  select(plot, rep, MAT, MDS1, MDS2)

#calculate tcrit, df = 7 for 95% CI 
t<- qt(0.975, df = 7)

#summarize wUF.MAT by calculating centroid of NMS scores and 95% CI for each MAT plot
centroid.wUF <- wUF.MAT %>%
  group_by(plot, MAT) %>%
  summarize(n = length(MDS1), 
    NMS1 = mean(MDS1), 
    NMS1.ci = t*(sd(MDS1)/sqrt(n)),
    NMS2 = mean(MDS2), 
    NMS2.ci = t*(sd(MDS2)/sqrt(n))
    )

#create custom legend title
leg.title <- expression(paste('MAT (', degree, 'C)'))

#ggplot2 code for ordination figure
p <- ggplot(centroid.wUF, aes(NMS1, NMS2)) +
  geom_errorbar(aes(ymax = NMS2 + NMS2.ci, ymin = NMS2 - NMS2.ci),
    width = 0.0, size = 0.22, colour = 'grey') +
  geom_errorbarh(aes(xmax = NMS1 + NMS1.ci, xmin = NMS1 - NMS1.ci),
    height = 0.0, size = 0.22, colour = 'grey') +
  geom_point(size = 3.8, aes(colour = MAT)) +
  theme_classic() +
    theme(axis.ticks = element_blank(), axis.text = element_blank()) +
    scale_colour_gradient(low = "black", high = "yellow", name = leg.title) +
    scale_x_continuous(limits = c(-0.28, 0.24)) +
  scale_y_continuous(limits = c(-0.15, 0.15)) +
  xlab('NMDS Axis 1') +
  ylab('NMDS Axis 2')

#save ordination figure as .pdf file
ggsave('wUF_NMDS.pdf', width = 5, height = 3.5)
#save ordination figure as .TIFF file
ggsave('wUF_NMDS.tiff', width = 5, height = 3.5, dpi = 400)

