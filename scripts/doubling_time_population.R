#!/usr/bin/env Rscript

##############################################################
#  script: doubling_time_population.R
#  author: Jing LI (GitHub ID: lj1221)
#  last edited: 2018.10.11
#  description:	plot doubling time of the populations evolving in the drugs across multiple time points
#  example: Rscript --vanilla doubling_time_population.R
##############################################################


library(plyr)
library(ggplot2)

script_path <- "/Users/litilab1/Downloads/LJ/Liti_Lab/PROJECTS/results/my_4way_data/GitHub_dir/four_parent_experimental_evolution_2017/scripts"
setwd(script_path)
data_path <- "./../phenotyping"

# doubling time of populations evolved in HU across all the time points (two-parent and four-parent populations)
data_HU_population <- read.table(paste0(data_path, "/doubling_time_population_HU.tsv"), header = TRUE, na.strings = "")
data_HU_population$doubling_time[data_HU_population$doubling_time == "NA"] <- NA
data_HU_population$doubling_time <- as.numeric(as.character(data_HU_population$doubling_time))
cdata_HU_population <- ddply (data_HU_population, c("condition_evolved", "condition_measurement",
                                                    "population", "cross", "replicate", "transfer"), summarise,
                              N = length(doubling_time),
                              mean = mean(doubling_time, na.rm = TRUE),
                              sd = sd(doubling_time, na.rm = TRUE)) 

cdata_HU_population$transfer = factor (cdata_HU_population$transfer, levels = c("P0","P1","P2","P3","P4",
                                                                                "P5","P6","P7","P8","P9",
                                                                                "P10","P11","P12","P13","P14",
                                                                                "P15", "P16", "P17"))
cdata_HU_population$population = factor (cdata_HU_population$population, levels = c("two-parent", "four-parent",
                                                                                    "WA", "NA", "WE", "SA"))


ggplot(subset(cdata_HU_population, (population == "two-parent" | population == "four-parent") & 
                (condition_measurement =="YNB+HU" & condition_evolved == "HU") & (transfer != "P17")),aes(x = transfer, y = mean, color = population))+
  geom_boxplot(lwd=0.5)+
  ggtitle("Hydroxyurea")+ylab("doubling time (h)")+
  xlab("transfer")+theme(axis.text.y  = element_text(size=20),axis.text.x  = element_text(size=20,angle=45))+
  theme(axis.title.x = element_text(size=20),axis.title.y = element_text(size=20))+
  theme(plot.title=element_text(face="bold", size=20))+
  theme(legend.title=element_blank(),legend.text = element_text(size=15))+
  theme(strip.text.x = element_text(size=20,face="bold"),strip.text.y = element_text(size=20, face="bold"))+
  ylim(0,8)+scale_color_brewer(palette="Set1")+theme_bw(base_size = 20)
ggsave("two-parent_four-parnt_HU.pdf", path = "./../results",
       scale = 1, width = 10, height = 5, dpi = 300)


# doubling time of populations evolved in RM across all the time points (two-parent and four-parent populations)
data_RM_population <- read.table(paste0(data_path, "/doubling_time_population_RM.tsv"), header = TRUE, na.strings = "")
data_RM_population$doubling_time[data_RM_population$doubling_time == "NA"] <- NA
data_RM_population$doubling_time <- as.numeric(as.character(data_RM_population$doubling_time))
cdata_RM_population <- ddply (data_RM_population, c("condition_evolved", "condition_measurement",
                                                    "population", "cross", "replicate", "transfer"), summarise,
                              N = length(doubling_time),
                              mean = mean(doubling_time, na.rm = TRUE),
                              sd = sd(doubling_time, na.rm = TRUE)) 

cdata_RM_population$transfer = factor (cdata_RM_population$transfer, levels = c("P0","P1","P2","P3","P4",
                                                                                "P5","P6","P7","P8","P9",
                                                                                "P10","P11","P12","P13","P14",
                                                                                "P15", "P16", "P17"))
cdata_RM_population$population = factor (cdata_RM_population$population, levels = c("two-parent", "four-parent",
                                                                                    "WA", "NA", "WE", "SA"))


ggplot(subset(cdata_RM_population, (population == "two-parent" | population == "four-parent") & 
                      (condition_measurement =="YNB+RM" & condition_evolved == "RM") & (transfer != "P17")),aes(x = transfer, y = mean, color = population))+
        geom_boxplot(lwd=0.5)+
        ggtitle("Rapamycin")+ylab("doubling time (h)")+
        xlab("transfer")+theme(axis.text.y  = element_text(size=20),axis.text.x  = element_text(size=20,angle=45))+
        theme(axis.title.x = element_text(size=20),axis.title.y = element_text(size=20))+
        theme(plot.title=element_text(face="bold", size=20))+
        theme(legend.title=element_blank(),legend.text = element_text(size=15))+
        theme(strip.text.x = element_text(size=20,face="bold"),strip.text.y = element_text(size=20, face="bold"))+
        ylim(0,8)+scale_color_brewer(palette="Set1")+theme_bw(base_size = 20)
ggsave("two-parent_four-parnt_RM.pdf", path = "./../results",
       scale = 1, width = 10, height = 5, dpi = 300)

