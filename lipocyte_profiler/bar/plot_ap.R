library(data.table)
library(tidyverse)
#library(treemap)
setwd("~/medusa/papers/TWAS/SNX10/lipocyte_profiler/bar")
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
dat<-as.data.frame(fread("ap_day14_sc_both_sexes.txt"),stringsAsFactors=FALSE)
dat<-dat[dat$q<0.1,]     
dat$color<-dat$feature_color
dat$color[grepl("Mito",dat$feature) & grepl("AGP",dat$feature)]<-"Combined Features"
dat$color[grepl("DNA",dat$feature) & grepl("BODIPY",dat$feature)]<-"Combined Features"
dat$color[grepl("AGP",dat$feature) & grepl("BODIPY",dat$feature)]<-"Combined Features"

#Munge data to make treemap
plot_dat<-as.data.frame(table(dat$color))
colnames(plot_dat)<-c("Feature","Num")

#Rename labels for plotting
plot_dat$Label<-NA
plot_dat$Label[plot_dat$Feature=="Combined Features"]<-"Combined features"
plot_dat$Label[plot_dat$Feature=="AGP_feature"]<-"Actin features"
plot_dat$Label[plot_dat$Feature=="bodipy_feature"]<-"Lipid features"
plot_dat$Label[plot_dat$Feature=="DNA_feature"]<-"DNA features"
plot_dat$Label[plot_dat$Feature=="Mito_feature"]<-"Mitochondrial features"
plot_dat$Label[plot_dat$Feature=="AGP_feature"]<-"Actin features"
plot_dat$Label[plot_dat$Feature=="other_feature"]<-"Other features"
plot_dat$Label=paste(plot_dat$Label,",\n N=",plot_dat$Num,sep='')


#Plot treemap
#pdf("SNX10_rs1534696_treemap.pdf",width=5.5,height=2.5)
#treemap(plot_dat,index="Label",vSize="Num",type="index",
#        palette=c("#C1C4CB",
#                  rgb(pony_colors[2,1:3]),
#                  rgb(pony_colors[7,1:3]),
#                  rgb(pony_colors[16,1:3]),
#                  rgb(pony_colors[11,1:3]),
#                  "#FFFFFF"),
#        fontcolor.labels="black",
#        align.labels=c("left","top"))
#dev.off()

#Plot poechart
pie(plot_dat$Num,labels=plot_dat$Num,col=c(rgb(pony_colors[2,1:3]),rgb(pony_colors[16,1:3]),"#C1C4CB",rgb(pony_colors[7,1:3]),rgb(pony_colors[11,1:3]),"#FFFFFF"))

#Plot bargraph
pdf("SNX10_rs1534696_AP_barplot.pdf",width=6.25,height=3.5)
ggplot()+
  geom_bar(data=plot_dat,aes(sapply(strsplit(Label,' '),'[[',1),y=Num,fill=Label),color="black",stat="identity")+
  scale_fill_manual(values=c(rgb(pony_colors[2,1:3]),"#C1C4CB",rgb(pony_colors[7,1:3]),rgb(pony_colors[16,1:3]),rgb(pony_colors[11,1:3]),"#FFFFFF"))+
  theme_minimal()+
  scale_y_continuous(name="Number of Features")+
  scale_x_discrete(name="Feature Class")+
  coord_flip()+
  theme(legend.position="none",
        axis.title=element_text(size=14),
        axis.text=element_text(size=14))
dev.off()

##################### Female ########################

#Read in significant AP features
dat<-as.data.frame(fread("ap_day14_sc_female.txt"),stringsAsFactors=FALSE)
dat<-dat[dat$q<0.1,]     
dat$color<-dat$feature_color
dat$color[grepl("Mito",dat$feature) & grepl("AGP",dat$feature)]<-"Combined Features"
dat$color[grepl("DNA",dat$feature) & grepl("BODIPY",dat$feature)]<-"Combined Features"
dat$color[grepl("AGP",dat$feature) & grepl("BODIPY",dat$feature)]<-"Combined Features"

#Munge data to make treemap
plot_dat<-as.data.frame(table(dat$color))
colnames(plot_dat)<-c("Feature","Num")

#Rename labels for plotting
plot_dat$Label<-NA
plot_dat$Label[plot_dat$Feature=="Combined Features"]<-"Combined features"
plot_dat$Label[plot_dat$Feature=="AGP_feature"]<-"Actin features"
plot_dat$Label[plot_dat$Feature=="bodipy_feature"]<-"Lipid features"
plot_dat$Label[plot_dat$Feature=="DNA_feature"]<-"DNA features"
plot_dat$Label[plot_dat$Feature=="Mito_feature"]<-"Mitochondrial features"
plot_dat$Label[plot_dat$Feature=="AGP_feature"]<-"Actin features"
plot_dat$Label[plot_dat$Feature=="other_feature"]<-"Other features"
plot_dat$Label=paste(plot_dat$Label,",\n N=",plot_dat$Num,sep='')


#Plot treemap
#pdf("SNX10_rs1534696_treemap.pdf",width=5.5,height=2.5)
#treemap(plot_dat,index="Label",vSize="Num",type="index",
#        palette=c("#C1C4CB",
#                  rgb(pony_colors[2,1:3]),
#                  rgb(pony_colors[7,1:3]),
#                  rgb(pony_colors[16,1:3]),
#                  rgb(pony_colors[11,1:3]),
#                  "#FFFFFF"),
#        fontcolor.labels="black",
#        align.labels=c("left","top"))
#dev.off()

#Plot poechart
pie(plot_dat$Num,labels=plot_dat$Num,col=c(rgb(pony_colors[2,1:3]),rgb(pony_colors[16,1:3]),"#C1C4CB",rgb(pony_colors[7,1:3]),rgb(pony_colors[11,1:3]),"#FFFFFF"))

#Plot bargraph
pdf("SNX10_rs1534696_AP_barplot_female.pdf",width=6.25,height=3.5)
ggplot()+
  geom_bar(data=plot_dat,aes(sapply(strsplit(Label,' '),'[[',1),y=Num,fill=Label),color="black",stat="identity")+
  scale_fill_manual(values=c(rgb(pony_colors[2,1:3]),"#C1C4CB",rgb(pony_colors[7,1:3]),rgb(pony_colors[16,1:3]),rgb(pony_colors[11,1:3]),"#FFFFFF"))+
  theme_minimal()+
  scale_y_continuous(name="Number of Features")+
  scale_x_discrete(name="Feature Class")+
  coord_flip()+
  theme(legend.position="none",
        axis.title=element_text(size=14),
        axis.text=element_text(size=14))
dev.off()
