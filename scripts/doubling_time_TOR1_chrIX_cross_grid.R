#!/usr/bin/env Rscript

##############################################################
#  script: doubling_time_isolates.R
#  author: Jing LI (GitHub ID: lj1221)
#  last edited: 2018.10.11
#  description:	plot doubling time of an array of diploids that all possible genotypes (TOR1 mutation/wildtype, chrIX one/two copies) were combined 
#  example: Rscript --vanilla doubling_time_TOR1_chrIX_cross_grid.R
##############################################################

library(plyr)
library(ggplot2)

script_path <- "/Users/litilab1/Downloads/LJ/Liti_Lab/PROJECTS/results/my_4way_data/GitHub_dir/four_parent_experimental_evolution_2017/scripts"
setwd(script_path)
data_path <- "./../phenotyping"

data_cross_grid <- read.table(paste0(data_path, "/doubling_time_cross_grid.tsv"), header = TRUE, na.strings = "")
data_cross_grid$doubling_time[data_cross_grid$doubling_time == "NA"] <- NA
data_cross_grid$doubling_time <- as.numeric(as.character(data_cross_grid$doubling_time))
data_cross_grid$crossed_strain<-factor(data_cross_grid$crossed_strain, levels=c("M1","M2","M3","M4",
                                                                      "M5","M6","M7","M8",
                                                                      "M9","M10","M11","M12",
                                                                      "M13","M14","M15","M16"))

ggplot(data_cross_grid,aes(x=crossed_strain,y=doubling_time,color=TOR1,shape=chrIX))+geom_boxplot(show.legend = F)+geom_point(size=2)+
  ggtitle("")+ylab("doubling time (hours)")+
  xlab("crossed strains")+theme(axis.text.y  = element_text(size=20),axis.text.x  = element_text(size=20,angle=45))+
  theme(axis.title.x = element_text(size=20),axis.title.y = element_text(size=20))+
  theme(plot.title=element_text(face="bold", size=20))+
  theme(legend.title=element_blank(),legend.text = element_text(size=15))+
  theme(strip.text.x = element_text(size=12,face="bold"),strip.text.y = element_text(size=12, face="bold"))+
  scale_color_brewer(palette="Set1")+theme_bw(base_size=20)
ggsave("TOR1_chrIX_cross_grid.pdf", path = "./../results",
       scale = 1, width = 10, height = 5, dpi = 300)


