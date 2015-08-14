# Paul Selmants
# January 2, 2014, modified February 3, 2015
## MAT bacteria OTU rarefaction curves ##

#load required packages
require(tidyr)
require(dplyr)
require(ggplot2)
require(RColorBrewer)

#read rarefaction data into R
data <- read.delim('observed_species.txt') 
#convert df from 'wide' to 'long' format using tidyr
#add MAT column using dplyr
#average OTUs across iterations within rep and sequencing depth 
rarefy <- data %>%
	gather(sampleID, OTUs, s934.8:s1024.3) %>%
	separate(sampleID, into = c('plot', 'rep'), sep = '\\.') %>%
	arrange(plot, rep) %>%
	mutate(MAT = c(rep(13.8, 880), rep(13.0, 660), rep(16.7, 880), rep('16.1b', 880), 
		rep(18.2, 880), rep(17.3, 880),rep('16.1a', 880), rep(15.5, 880), rep(15.1, 770))) %>% 
	group_by(MAT, rep, sequences.per.sample) %>%
	summarise(iterations = length(OTUs), 
			OTU = mean(OTUs)) 

#calculate mean and se of OTUs by rep in each MAT plot 
rarefy.summary <- rarefy %>%
	group_by(MAT, sequences.per.sample) %>%
	summarise(N = length(OTU), 
		OTUs = mean(OTU), 
		sd   = sd(OTU), 
		se   = sd/sqrt(N))

#Custom legend title
leg.title <- expression(paste('MAT (', degree, 'C)'))

#Make graph of rarefaction curves
p1 <- ggplot(rarefy.summary, aes(sequences, OTUs, color = MAT, group = MAT)) + 
	geom_point(size = 3) + 
	theme_classic() +
	scale_y_continuous(limits = c(0,800)) +
	geom_smooth(method = loess, se = FALSE, size = 0.3) +
	geom_errorbar(aes(ymax = OTUs + se, ymin = OTUs - se), width = 0, size = 0.2) +
	scale_colour_manual(values = rev(brewer.pal(9, 'Spectral')), name = leg.title) +
	ylab('observed OTUs') +
	xlab('Sequences per sample') 

#save figure as .pdf
ggsave('raref_lines2.pdf', width = 7, height = 4)

