#!/usr/bin/env Rscript

##############################################################
#  script: doubling_time_isolates.R
#  author: Jing LI (GitHub ID: lj1221)
#  last edited: 2018.10.11
#  description:	plot doubling time of the TOR1/2 variants on various nitrogen limited conditions
#  example: Rscript --vanilla doubling_time_TOR_variants_nitrogen_limit.R
##############################################################

library(plyr)
library(ggplot2)
library(RColorBrewer)

script_path <- "/Users/litilab1/Downloads/LJ/Liti_Lab/PROJECTS/results/my_4way_data/GitHub_dir/four_parent_experimental_evolution_2017/scripts"
setwd(script_path)
data_path <- "./../phenotyping"

data_TOR <- read.table(paste0(data_path, "/doubling_time_TOR_variants_in_nitrogen_limited_conditions.tsv"), header = TRUE, na.strings = "")
data_TOR$Strain = factor(data_TOR$Strain, levels=c("WE/SA_TOR2Del_YGL2915", "WE_TOR2Del/SA_YGL2914",
                                                   "WE/SA_TOR1Del_YGL2498","WE_TOR1Del/SA_YGL2497", 
                                                   "SA/SA_CC454","WE/WE_CC411"))
myPalette <- colorRampPalette(rev(brewer.pal(100, "RdBu")))
myPalette <- colorRampPalette(c("red","white","blue"))
ggplot(subset(data_TOR,(Strain!="SA/WE_CC444") & 
                (Condition!="YNB_complete+RM" & Condition!="YNB_complete+D2O")),
       aes(x=Condition,y=Strain)) + geom_tile(aes(fill=doubling_time_compared_with_SA.WE)) + 
  coord_equal()+
  scale_fill_gradientn(colours= myPalette(100),na.value = "transparent")+
  ggtitle("")+ylab("strain")+
  xlab("condition")+theme(axis.text.y  = element_text(size=15),axis.text.x  = element_text(size=10,angle=90))+
  theme(axis.title.x = element_text(size=20),axis.title.y = element_text(size=20))+
  theme(plot.title=element_text(face="bold", size=20))+
  theme(legend.title=element_blank(),legend.text = element_text(size=15))+
  theme(strip.text.x = element_text(size=12,face="bold"),strip.text.y = element_text(size=12, face="bold"))+
  theme(axis.text.x = element_text(margin = margin(0, unit = "cm")))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank())
ggsave("TOR1-2_variants_nitrogen_limited_condition.pdf", path = "./../results",
       scale = 1, width = 10, height = 10, dpi = 300)

