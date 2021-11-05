library(data.table)
library(tidyverse)
#library(treemap)
setwd("~/medusa/papers/TWAS/lipocyte_profiler/bar")
pony_colors<-fread("~/medusa/papers/TWAS/pony_palette.txt")

########## For color manipulation of graph ###############
darken <- function(color, factor=1.1){
  col <- col2rgb(color)
  col <- col/factor
  col <- rgb(t(col), maxColorValue=255)
  col
}
lighten <- function(color, factor=0.8){
  col <- col2rgb(color)
  col <- col/factor
  for (i in 1:length(col)) { if ( col[i] > 255 ) { col[i] = 255 } }
  col <- rgb(t(col), maxColorValue=255)
  col
}
#########################################################

##################### Both sexes ########################

#Read in significant AP features
all<-as.data.frame(fread("../scatter/rs1534696_allfeatures_day14_sc_bothsexes_nonredundant.txt"),stringsAsFactors=FALSE)
sig<-all[all$q<0.1 & all$SNP.pvalue<0.05,]
sig$color<-NA
sig$color[grepl("Mito",sig$feature)]<-"Mito_feature"
sig$color[grepl("DNA",sig$feature)]<-"DNA_feature"
sig$color[grepl("AGP",sig$feature)]<-"AGP_feature"
sig$color[grepl("BODIPY",sig$feature)]<-"bodipy_feature"
sig$color[grepl("Mito",sig$feature) & grepl("AGP",sig$feature)]<-"Combined Features"
sig$color[grepl("DNA",sig$feature) & grepl("BODIPY",sig$feature)]<-"Combined Features"
sig$color[grepl("AGP",sig$feature) & grepl("BODIPY",sig$feature)]<-"Combined Features"
sig$color[grepl("Mito",sig$feature) & grepl("BODIPY",sig$feature)]<-"Combined Features"
sig$color[grepl("Mito",sig$feature) & grepl("DNA",sig$feature)]<-"Combined Features"
sig$color[grepl("DNA",sig$feature) & grepl("AGP",sig$feature)]<-"Combined Features"
sig$color[is.na(sig$color)]<-"other_feature"

#Read in all features, color by category
all<-as.data.frame(fread("../scatter/rs1534696_allfeatures_day14_sc_bothsexes_nonredundant.txt"))
all$color<-NA
all$color[grepl("Mito",all$feature)]<-"Mito_feature"
all$color[grepl("DNA",all$feature)]<-"DNA_feature"
all$color[grepl("AGP",all$feature)]<-"AGP_feature"
all$color[grepl("BODIPY",all$feature)]<-"bodipy_feature"
all$color[grepl("Mito",all$feature) & grepl("AGP",all$feature)]<-"Combined Features"
all$color[grepl("DNA",all$feature) & grepl("BODIPY",all$feature)]<-"Combined Features"
all$color[grepl("AGP",all$feature) & grepl("BODIPY",all$feature)]<-"Combined Features"
all$color[grepl("Mito",all$feature) & grepl("BODIPY",all$feature)]<-"Combined Features"
all$color[grepl("Mito",all$feature) & grepl("DNA",all$feature)]<-"Combined Features"
all$color[grepl("DNA",all$feature) & grepl("AGP",all$feature)]<-"Combined Features"
all$color[is.na(all$color)]<-"other_feature"

#Evaluate signifiance as: sig features in that category/all features in that category
plot_dat<-as.data.frame(t(table(sig$color)/table(all$color)))
plot_dat<-plot_dat[,2:3]
colnames(plot_dat)<-c("Feature","Proportion significant")

#Rename labels for plotting
plot_dat$Label<-NA
plot_dat$Label[plot_dat$Feature=="Combined Features"]<-"Combined features"
plot_dat$Label[plot_dat$Feature=="AGP_feature"]<-"Actin features"
plot_dat$Label[plot_dat$Feature=="bodipy_feature"]<-"Lipid features"
plot_dat$Label[plot_dat$Feature=="DNA_feature"]<-"DNA features"
plot_dat$Label[plot_dat$Feature=="Mito_feature"]<-"Mitochondrial features"
plot_dat$Label[plot_dat$Feature=="AGP_feature"]<-"Actin features"
plot_dat$Label[plot_dat$Feature=="other_feature"]<-"Other features"
plot_dat$Label=paste(plot_dat$Label,",\n N=",paste(table(sig$color),table(all$color),sep='/'),sep='')

#Plot bargraph: proportion of sig features in each category, colored by category
pdf("SNX10_rs1534696_AP_barplot_props.pdf",width=6.25,height=3.5)
ggplot()+
  geom_bar(data=plot_dat,aes(sapply(strsplit(Label,' '),'[[',1),y=`Proportion significant`,fill=Label),color="black",stat="identity")+
  scale_fill_manual(values=c(rgb(pony_colors[2,1:3]),"#C1C4CB",rgb(pony_colors[7,1:3]),rgb(pony_colors[16,1:3]),rgb(pony_colors[11,1:3]),"#FFFFFF"))+
  theme_minimal()+
  scale_x_discrete(name="Feature Class")+
  coord_flip()+
  theme(legend.position="none",
        axis.title=element_text(size=14),
        axis.text=element_text(size=14))
dev.off()
