# Paul Selmants
# March 26, 2015
## MAT bacteria phylum relative abundance ## 

#load required packages
require(tidyr)
require(dplyr)
require(ggplot2)
require(RColorBrewer)
#read phylum relative abundance data into R
data <- read.delim('Site_otu_table_sorted_L2.txt') 
#convert df from 'wide' to 'long' format using tidyr
#add MAT column using dplyr
phyla <- data %>%
	gather(plot, relabun, h1468:w1274) %>% 
	separate(Taxon, into = c('kingdom', 'phylum'), sep = '\\;p__') %>%
	mutate(relabun.pct = relabun*100,
		MAT = rep(c(13.8, 13.0, 18.2, 17.3, 16.7, '16.1b', '16.1a', 15.5, 15.1), each = 36))
#sum phyla with low relabun into 'Other' category
low.phyla <- phyla %>%
	filter(!grepl('Acidobacteria|Actinobacteria|Chloroflexi|Planctomycetes|Proteobacteria|WPS-2',
		phylum)) %>%
	group_by(MAT, plot) %>%
	summarise(relabun.pct = sum(relabun.pct)) %>%
	mutate(phylum = 'Other')
#filter out all but six most abundant phyla 
high.phyla <- phyla %>%
	filter(grepl('Proteobacteria|Acidobacteria|Actinobacteria|Chloroflexi|WPS-2|Planctomycetes',
		phylum)) %>%
	select(MAT, plot, relabun.pct, phylum)
#merge into one df and arrange by MAT
phylumfig <- rbind(high.phyla, low.phyla) %>%
	arrange(MAT, factor(phylum, levels = c('Proteobacteria', 'Acidobacteria', 'Actinobacteria',
		 'Chloroflexi' , 'WPS-2', 'Planctomycetes', 'Other')))
#colorblind-friendly palette for fig
cbPalette <- c("#0072B2", "#009E73", "#56B4E9", "#999999", "#F0E442","#D55E00", "#CC79A7")
#stacked bar chart figure 
p <- ggplot(phylumfig, aes(MAT, relabun.pct, fill = factor(phylum))) + 
	geom_bar(stat="identity") +
	theme_classic() +
	scale_y_continuous(expand = c(0,0)) +
	xlab(expression(paste('Mean annual temperature (', degree, 'C)'))) +
	theme(axis.title.x = element_text(vjust = 0.05)) + #axis title further from x-axis
	ylab("Relative abundance (%)") +
	scale_fill_manual(values = cbPalette, name='Phylum',
		breaks = c('Proteobacteria', 'Acidobacteria', 'Actinobacteria', 'Chloroflexi', 
		'WPS-2', 'Planctomycetes', 'Other'))

#save stacked bar chart as pdf
ggsave("MATphyla.pdf", width = 7, height = 3.5)
