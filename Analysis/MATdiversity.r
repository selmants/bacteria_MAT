#Paul Selmants
#March 31, 2015
##MAT bacterial diversity indices##

#load required packages
library(tidyr)
library(ggplot2)
library(RColorBrewer)
library(dplyr)

#read data into R
obs <- read.delim('observed_species.txt')
chao1 <- read.delim('chao1.txt')
pd <- read.delim('PD_whole_tree.txt')
#Use tidyr to convert OTU richness df from 'wide' to 'long' format
#Use dplyr to filter columns by 2200 sequencing depth and add MAT column
observed <- obs %>%
	gather(sampleID, OTUs, s934.8:s1024.3) %>%
	separate(sampleID, into = c('plot', 'rep'), sep = '\\.') %>%
	arrange(plot, rep) %>%
	mutate(MAT = c(rep(13.8, 880), rep(13.0, 660), rep(16.7, 880), rep(16.1, 880), 
		rep(18.2, 880), rep(17.3, 880),rep(16.1, 880), rep(15.5, 880), rep(15.1, 770))) %>%
	filter(sequences.per.sample == 2200) %>%
	group_by(plot, MAT, rep) %>%
	summarise(iterations = length(OTUs), 
			OTU.rich = mean(OTUs))
#Use tidyr to convert chao1 estimated OTU richness df from 'wide' to 'long' format
#Use dplyr to filter columns by 2200 sequencing depth and add MAT column
chao1.rich <- chao1 %>%
	gather(sampleID, chao1, s934.8:s1024.3) %>%
	separate(sampleID, into = c('plot', 'rep'), sep = '\\.') %>%
	arrange(plot, rep) %>%
	mutate(MAT = c(rep(13.8, 880), rep(13.0, 660), rep(16.7, 880), rep(16.1, 880), 
		rep(18.2, 880), rep(17.3, 880),rep(16.1, 880), rep(15.5, 880), rep(15.1, 770))) %>%
	filter(sequences.per.sample == 2200) %>%
	group_by(plot, MAT, rep) %>%
	summarise(iterations = length(chao1), 
			chao1.rich = mean(chao1))
#join observed and estimated (chao1) OTU richness estimates
richness <- full_join(observed, chao1.rich)
#Use tidyr to convert Faith's phylogenetic diversity df from 'wide' to 'long' format
#Use dplyr to filter columns by 2200 sequencing depth and add MAT column
phyl.div <- pd %>% 
	gather(sampleID, pd, s934.8:s1024.3) %>%
	separate(sampleID, into = c('plot', 'rep'), sep = '\\.') %>%
	arrange(plot, rep) %>%
	mutate(MAT = c(rep(13.8, 880), rep(13.0, 660), rep(16.7, 880), rep(16.1, 880), 
		rep(18.2, 880), rep(17.3, 880),rep(16.1, 880), rep(15.5, 880), rep(15.1, 770))) %>%
	filter(sequences.per.sample == 2200) %>%
	group_by(plot, MAT, rep) %>%
	summarise(iterations = length(pd), 
			pd = mean(pd))
#join phylogenetic diversity with richness estimates 
diversity <- full_join(richness, phyl.div) %>%
	select(-iterations) %>%
	arrange(MAT, plot)
#calculate tcrit value for df = 7
t<- qt(0.975, df = 7)
#summarize diversity estimates (mean and 95% CI)
div.summary <- diversity %>%
	group_by(plot, MAT) %>%
	summarise(N = length(OTU.rich), 
		obs = mean(OTU.rich),
		obs.ci = t*(sd(OTU.rich)/sqrt(N)), 
		chao1 = mean(chao1.rich),
		chao1.ci = t*(sd(chao1.rich)/sqrt(N)),
		FaithPD = mean(pd), 
		FaithPD.ci = t*(sd(pd)/sqrt(N)))

#Linear regression of Faith's PD as a function of MAT
pd.mod <- lm(FaithPD ~ MAT, data = div.summary)
#Linear regression of estimated OTU richness (Chao1) as a function of MAT
chao.mod <- lm(chao1 ~ MAT, data = div.summary)
#Linear regression of observed OTU richness as a function of MAT
obs.mod <- lm(obs ~ MAT, data = div.summary)


#custom x-axis label for figures
my.xlab = expression(paste('Mean annual temperature (', degree, 'C)'))

#use ggplot2 to make fig of pd vs. MAT
pd.fig <- ggplot(div.summary, aes(MAT, FaithPD)) + 
	theme_classic() +
	geom_point(size = 2.8, colour = 'blue') + 
	geom_errorbar(aes(ymax = FaithPD + FaithPD.ci, ymin = FaithPD - FaithPD.ci), 
		width = 0.04, size = 0.2, colour = 'blue') +
	scale_y_continuous(limits = c(0,80), breaks = c(0,20,40,60,80)) +
	xlab(my.xlab) +
	ylab("Phylogenetic diversity") +
	theme(axis.title.y = element_text(vjust = 1.25)) #axis title further from y-axis
#save figure as .pdf file
ggsave('FaithPD.pdf', width = 4, height = 3.5)

obs.fig <- ggplot(div.summary, aes(MAT, obs)) +
	theme_classic() + 
	geom_point(size = 2.8, colour = 'blue') +
	geom_errorbar(aes(ymax = obs + obs.ci, ymin = obs - obs.ci), 
		width = 0.04, size = 0.2, colour = 'blue') +
	scale_y_continuous(limits = c(0, 1000), breaks = c(0,250,500,750,1000)) +
	xlab(my.xlab) +
	ylab('Observed OTU richness') +
	theme(axis.title.y = element_text(vjust = 0.7))
#save figure as .pdf file
ggsave('ObsRich.pdf', width = 4, height = 3.5)

#save FaithPD and obsrich as two-panel figure in TIFF format
library(gridExtra)
tiff('PD_rich.tiff', width = 2500, height = 1050, res = 400)
grid.arrange(obs.fig, pd.fig, ncol = 2)
dev.off()

