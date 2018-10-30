#!/usr/bin/env Rscript

##############################################################
#  script: doubling_time_isolates.R
#  author: Jing LI (GitHub ID: lj1221)
#  last edited: 2018.10.11
#  description:	plot doubling time of the individuals evolving in the drugs at initial and end time points
#  example: Rscript --vanilla doubling_time_isolates.R
##############################################################


library(plyr)
library(ggplot2)

script_path <- "/Users/litilab1/Downloads/LJ/Liti_Lab/PROJECTS/results/my_4way_data/GitHub_dir/four_parent_experimental_evolution_2017/scripts"
setwd(script_path)
data_path <- "./../phenotyping"

#doubling time of isolates from four-parent populations and isogenic parental populations evolved in HU, the initial and end time points
data_HU_isolates <- read.table(paste0(data_path, "/doubling_time_isolates_HU.tsv"), header = TRUE, na.strings = "")
data_HU_isolates$doubling_time[data_HU_isolates$doubling_time == "NA"] <- NA
data_HU_isolates$doubling_time <- as.numeric(as.character(data_HU_isolates$doubling_time))
cdata_HU_isolates <- ddply (data_HU_isolates, c("condition_evolved", "condition_measurement",
                                                    "population", "isolate", 
                                                "genotyped_mutated_gene", "genotyped_mutated_site", "day"), summarise,
                              N = length(doubling_time),
                              mean = mean(doubling_time, na.rm = TRUE),
                              sd = sd(doubling_time, na.rm = TRUE))

cdata_HU_isolates_four_parent_T0 <- subset(cdata_HU_isolates,condition_evolved=="MO" & condition_measurement=="YNB+HU" & 
                                             (population=="F12_1_MO_1"|population=="F12_2_MO_1"))
cdata_HU_isolates_four_parent_end <- subset(cdata_HU_isolates,condition_evolved=="HU" & condition_measurement=="YNB+HU" & 
                                              (population=="F12_1_HU_1"|population=="F12_1_HU_2"|population=="F12_2_HU_2"|population=="F12_2_HU_3"))
cdata_HU_isolates_WA_T0 <- subset(cdata_HU_isolates,condition_evolved=="MO" & condition_measurement=="YNB+HU" & population=="WA")
cdata_HU_isolates_NA_T0 <- subset(cdata_HU_isolates,condition_evolved=="MO" & condition_measurement=="YNB+HU" & population=="NA")
cdata_HU_isolates_WE_T0 <- subset(cdata_HU_isolates,condition_evolved=="MO" & condition_measurement=="YNB+HU" & population=="WE")
cdata_HU_isolates_SA_T0 <- subset(cdata_HU_isolates,condition_evolved=="MO" & condition_measurement=="YNB+HU" & population=="SA")
cdata_HU_isolates_NA_end <- subset(cdata_HU_isolates,condition_evolved=="HU" & condition_measurement=="YNB+HU" & 
                                     (population=="NA_HU_1"|population=="NA_HU_2"))
cdata_HU_isolates_WE_end <- subset(cdata_HU_isolates,condition_evolved=="HU" & condition_measurement=="YNB+HU" & 
                                     (population=="WE_HU_1"|population=="WE_HU_2"))
cdata_HU_isolates_SA_end <- subset(cdata_HU_isolates,condition_evolved=="HU" & condition_measurement=="YNB+HU" & 
                                     (population=="SA_HU_1"|population=="SA_HU_2"))

ggplot(cdata_HU_isolates_four_parent_T0,aes(x="four-parent",y=mean, color="T0", alpha=0.1))+geom_boxplot(lwd=0.8)+
  geom_boxplot(data=cdata_HU_isolates_four_parent_end,aes(x="four-parent",y=mean, color="T14"),lwd=0.8)+
  geom_boxplot(data=cdata_HU_isolates_WA_T0,aes(x="WA",y=mean, color="T0"),lwd=0.8)+
  geom_boxplot(data=cdata_HU_isolates_NA_T0,aes(x="NA",y=mean, color="T0"),lwd=0.8)+
  geom_boxplot(data=cdata_HU_isolates_NA_end,aes(x="NA",y=mean, color="T14"),lwd=0.8)+
  geom_boxplot(data=cdata_HU_isolates_WE_T0,aes(x="WE",y=mean, color="T0"),lwd=0.8)+
  geom_boxplot(data=cdata_HU_isolates_WE_end,aes(x="WE",y=mean, color="T14"),lwd=0.8)+
  geom_boxplot(data=cdata_HU_isolates_SA_T0,aes(x="SA",y=mean, color="T0"),lwd=0.8)+
  geom_boxplot(data=cdata_HU_isolates_SA_end,aes(x="SA",y=mean, color="T14"),lwd=0.8)+
  scale_x_discrete(limits=c("four-parent","WA","NA","WE","SA"))+
  ylim(0,8)+
  scale_color_manual(values=c("#1B9E77", "#D95F02"))+
  ggtitle("isolates HU")+ylab("doubling time (h)")+xlab("")+
  theme(axis.text.x  = element_text(size=20),axis.text.y  = element_text(size=20))+
  theme(axis.title.x = element_text(size=20),axis.title.y = element_text(size=20))+
  theme(plot.title=element_text(face="bold", size=20))+
  theme(legend.title=element_blank(),legend.text = element_text(size=15))+
  theme(strip.text.x = element_text(size=20,face="bold"),strip.text.y = element_text(size=20, face="bold"))+
  theme_bw(base_size=20)
ggsave("isolates_HU.pdf", path = "./../results",
       scale = 1, width = 10, height = 5, dpi = 300)



#doubling time of isolates from four-parent populations and isogenic parental populations evolved in RM, the initial and end time points
data_RM_isolates <- read.table(paste0(data_path, "/doubling_time_isolates_RM.tsv"), header = TRUE, na.strings = "")
data_RM_isolates$doubling_time[data_RM_isolates$doubling_time == "NA"] <- NA
data_RM_isolates$doubling_time <- as.numeric(as.character(data_RM_isolates$doubling_time))
cdata_RM_isolates <- ddply (data_RM_isolates, c("condition_evolved", "condition_measurement",
                                                "population", "isolate", 
                                                "genotyped_mutated_gene", "genotyped_mutated_site", "day"), summarise,
                            N = length(doubling_time),
                            mean = mean(doubling_time, na.rm = TRUE),
                            sd = sd(doubling_time, na.rm = TRUE))

cdata_RM_isolates_four_parent_T0 <- subset(cdata_RM_isolates,condition_evolved=="MO" & condition_measurement=="YNB+RM" & 
                                             (population=="F12_1_MO_1"|population=="F12_2_MO_1"))
cdata_RM_isolates_four_parent_end <- subset(cdata_RM_isolates,condition_evolved=="RM" & condition_measurement=="YNB+RM" & 
                                              (population=="F12_1_RM_2"|population=="F12_1_RM_4"|population=="F12_2_RM_2"|population=="F12_2_RM_4"))
cdata_RM_isolates_WA_T0 <- subset(cdata_RM_isolates,condition_evolved=="MO" & condition_measurement=="YNB+RM" & population=="WA")
cdata_RM_isolates_NA_T0 <- subset(cdata_RM_isolates,condition_evolved=="MO" & condition_measurement=="YNB+RM" & population=="NA")
cdata_RM_isolates_WE_T0 <- subset(cdata_RM_isolates,condition_evolved=="MO" & condition_measurement=="YNB+RM" & population=="WE")
cdata_RM_isolates_SA_T0 <- subset(cdata_RM_isolates,condition_evolved=="MO" & condition_measurement=="YNB+RM" & population=="SA")
cdata_RM_isolates_WA_end <- subset(cdata_RM_isolates,condition_evolved=="RM" & condition_measurement=="YNB+RM" & 
                                     (population=="WA_RM_1"|population=="WA_RM_2"))
cdata_RM_isolates_NA_end <- subset(cdata_RM_isolates,condition_evolved=="RM" & condition_measurement=="YNB+RM" & 
                                     (population=="NA_RM_1"|population=="NA_RM_2"))
cdata_RM_isolates_WE_end <- subset(cdata_RM_isolates,condition_evolved=="RM" & condition_measurement=="YNB+RM" & 
                                     (population=="WE_RM_1"|population=="WE_RM_2"))
cdata_RM_isolates_SA_end <- subset(cdata_RM_isolates,condition_evolved=="RM" & condition_measurement=="YNB+RM" & 
                                     (population=="SA_RM_1"|population=="SA_RM_2"))

ggplot(cdata_RM_isolates_four_parent_T0,aes(x="four-parent",y=mean, color="T0", alpha=0.1))+geom_boxplot(lwd=0.8)+
  geom_boxplot(data=cdata_RM_isolates_four_parent_end,aes(x="four-parent",y=mean, color="T15"),lwd=0.8)+
  geom_boxplot(data=cdata_RM_isolates_WA_T0,aes(x="WA",y=mean, color="T0"),lwd=0.8)+
  geom_boxplot(data=cdata_RM_isolates_WA_end,aes(x="WA",y=mean, color="T15"),lwd=0.8)+
  geom_boxplot(data=cdata_RM_isolates_NA_T0,aes(x="NA",y=mean, color="T0"),lwd=0.8)+
  geom_boxplot(data=cdata_RM_isolates_NA_end,aes(x="NA",y=mean, color="T15"),lwd=0.8)+
  geom_boxplot(data=cdata_RM_isolates_WE_T0,aes(x="WE",y=mean, color="T0"),lwd=0.8)+
  geom_boxplot(data=subset(cdata_RM_isolates_WE_end,population=="WE_RM_1"),aes(x="WE",y=mean, color="T15"),lwd=0.8)+
  geom_boxplot(data=subset(cdata_RM_isolates_WE_end,population=="WE_RM_2"),aes(x="WE",y=mean, color="T8"),lwd=0.8)+
  geom_boxplot(data=cdata_RM_isolates_SA_T0,aes(x="SA",y=mean, color="T0"),lwd=0.8)+
  geom_boxplot(data=cdata_RM_isolates_SA_end,aes(x="SA",y=mean, color="T15"),lwd=0.8)+
  scale_x_discrete(limits=c("four-parent","WA","NA","WE","SA"))+ylim(0,8)+
  scale_color_manual(values=c("#1B9E77", "#D95F02","#E7298A"))+
  ggtitle("isolates RM")+ylab("doubling time (h)")+xlab("")+
  theme(axis.text.x  = element_text(size=20),axis.text.y  = element_text(size=20))+
  theme(axis.title.x = element_text(size=20),axis.title.y = element_text(size=20))+
  theme(plot.title=element_text(face="bold", size=20))+
  theme(legend.title=element_blank(),legend.text = element_text(size=15))+
  theme(strip.text.x = element_text(size=20,face="bold"),strip.text.y = element_text(size=20, face="bold"))+
  theme_bw(base_size=20)
ggsave("isolates_RM.pdf", path = "./../results",
       scale = 1, width = 10, height = 5, dpi = 300)