#Pielou's J (evenness) for MAT bacteria data
#Paul Semants
#May 15, 2015

library(vegan)
library(dplyr)
library(tidyr)

data <- read.delim('Site_otu_table_sorted_L6.txt', row.names = 1)

tdata <- t(data)

H <-diversity(tdata)
J <- H/log(specnumber(tdata))

Pielou <- as.data.frame(J) 

PielouJ <- setNames(cbind(rownames(Pielou), Pielou, row.names = NULL), 
         c("plotID", "Jscore")) 

PielouJ <- PielouJ %>%
	mutate(MAT = c(13.8, 13.0, 18.2, 17.3, 16.7, 16.1, 16.1, 15.5, 15.1))
	
