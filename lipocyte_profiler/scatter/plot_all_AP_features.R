#!/usr/bin/Rscript
library(data.table)
library(tidyverse)
library(gridExtra)
setwd("~/medusa/papers/TWAS/lipocyte_profiler/scatter")
pony_colors<-fread("~/medusa/papers/TWAS/pony_palette")

######### Author: Grace Hansen #########
#This script plots adipocyte profiler data along a time course and from different cell types

############################ For color manipulation ############################
darken <- function(color, factor=1.2){
  col <- col2rgb(color)
  col <- col/factor
  col <- rgb(t(col), maxColorValue=255)
  col
}
################################################################################
celltypes<-c("sc","vc")
timepoints<-c("day0","day3","day8","day14")

################# Both sexes ####################3
AP<-matrix(nrow=0,ncol=9)
for (ct in celltypes) {
  for (tp in timepoints) {
    dat<-as.matrix(fread(paste("rs1534696_allfeatures_",tp,"_",ct,"_bothsexes_nonredundant.txt",sep='')))
    head(dat)
    dat<-cbind(dat,rep(ct,nrow(dat)))
    dat<-cbind(dat,rep(tp,nrow(dat)))
    AP<-rbind(AP,dat)
  }
}
AP<-as.data.frame(AP,stringsAsFactors=FALSE)
colnames(AP)[8:9]<-c("celltype","tp")
AP$`SNP.pvalue`<-as.numeric(as.character(AP$`SNP.pvalue`))
AP$`t-test`<-as.numeric(as.character(AP$`t-test`))
AP$q<-as.numeric(as.character(AP$q))
AP$tp<-factor(AP$tp,levels=c("day0","day3","day8","day14"))

#Insert category labels
AP$category<-rep("other/combined",nrow(AP))
AP$category[grepl("AGP",AP$features) & !(grepl("BODIPY",AP$features)) & !(grepl("Mito",AP$features)) & !(grepl("DNA",AP$features))]<-"Actin"
AP$category[!(grepl("AGP",AP$features)) & !(grepl("BODIPY",AP$features)) & !(grepl("Mito",AP$features)) & grepl("DNA",AP$features)]<-"DNA"
AP$category[!(grepl("AGP",AP$features)) & grepl("BODIPY",AP$features) & !(grepl("Mito",AP$features)) & !(grepl("DNA",AP$features))]<-"Intracellular lipids"
AP$category[!(grepl("AGP",AP$features)) & !(grepl("BODIPY",AP$features)) & grepl("Mito",AP$features) & !(grepl("DNA",AP$features))]<-"Mitochondria"
AP<-AP[!(is.na(AP$SNP.pvalue)),]

#Volcano plots
S<-ggplot()+
  geom_hline(yintercept=1.855,linetype="dashed",color="gray80")+
  geom_point(data=AP[AP$celltype=="sc" & (AP$q > 0.1 |  AP$`SNP.pvalue`>0.05),],aes(x=`t-test`,y=-log10(`SNP.pvalue`),color=category),size=0.5,alpha=0.25)+
  geom_point(data=AP[AP$celltype=="sc" & AP$q < 0.1 & AP$`SNP.pvalue`<0.05,],aes(x=`t-test`,y=-log10(`SNP.pvalue`),color=category),size=1.5)+
  facet_wrap(vars(tp),nrow=1)+
  ggtitle("Subcutaneous Adipocytes")+
  theme_minimal()+
  scale_x_continuous(name="t-statistic",limits=c(-4,4))+
  scale_y_continuous(name="-log10 p-value",limits=c(0,5))+
  scale_color_manual(values=c(rgb(pony_colors[2,1:3]),rgb(pony_colors[7,1:3]),darken(rgb(pony_colors[16,1:3]),1.1),rgb(pony_colors[11,1:3]),"gray80"))+
  theme(axis.line.x = element_blank(),
        axis.text.x=element_text(size=8))+
  labs(color="Feature Class")

V<-ggplot()+
  geom_hline(yintercept=1.855,linetype="dashed",color="gray80")+
  geom_point(data=AP[AP$celltype=="vc" & (AP$q > 0.1 |  AP$`SNP.pvalue`>0.05),],aes(x=`t-test`,y=-log10(`SNP.pvalue`),color=category),size=0.5,alpha=0.25)+
  geom_point(data=AP[AP$celltype=="vc" & AP$q < 0.1 & AP$`SNP.pvalue`<0.05,],aes(x=`t-test`,y=-log10(`SNP.pvalue`),color=category),size=1.5)+
  facet_wrap(vars(tp),nrow=1)+
  ggtitle("Visceral Adipocytes")+
  theme_minimal()+
  scale_x_continuous(name="t-statistic",limits=c(-4,4))+
  scale_y_continuous(name="-log10 p-value",limits=c(0,5))+
  scale_color_manual(values=c(rgb(pony_colors[2,1:3]),rgb(pony_colors[7,1:3]),darken(rgb(pony_colors[16,1:3]),1.1),rgb(pony_colors[11,1:3]),"gray80"))+
  theme(axis.line.x = element_blank(),
        axis.text.x=element_text(size=8))+
  labs(color="Feature Class")

pdf("AP_features_timecourse.pdf",width=8,height=3)
grid.arrange(V,S,nrow=1)
dev.off()
