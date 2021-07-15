#!/usr/bin/Rscript
library(tidyverse)
library(data.table)
library(gridExtra)

trait<-"WHR"
setwd(paste("/home/grace/midway/ldsc_seg",sep=""))
pony_colors<-fread("~/medusa/papers/TWAS/pony_palette")

######### Author: Grace Hansen #########
#This script plots adipocyte/preadipocyte cell proportion estimates from xCell for GTEx tissues. 

#Load and plot ldsc-seg heritability estimates
dat<-fread(paste("output/",trait,"/GTEx_exp/all_tissues.results",sep=''))
groups<-fread(paste("output/",trait,"/GTEx_exp/tissues_groups",sep=''))
dat<-merge(groups,dat,by="Tissue")
dat<-dat %>% arrange(Group,Enrichment)
dat$Tissue<-factor(dat$Tissue,levels=unique(dat$Tissue))

set.seed(29)
pony_palette=slice(pony_colors,sample(1:length(unique(dat$Group))))

#Load xCell estimates, generate fat content
prep_columns<-function(dat) {
  CD4<-colMeans(dat[c("CD4+ memory T-cells","CD4+ naive T-cells","CD4+ T-cells","CD4+ Tcm","CD4+ Tem"),])
  CD8<-colMeans(dat[c("CD8+ naive T-cells","CD8+ T-cells","CD8+ Tcm","CD8+ Tem"),])
  DC<-colMeans(dat[c("aDC","cDC","DC","iDC","pDC"),])
  Endothelial<-colMeans(dat[c("Endothelial cells","ly Endothelial cells","mv Endothelial cells"),])
  Macrophage<-colMeans(dat[c("Macrophages","Macrophages M1","Macrophages M2"),])
  Bcell<-colMeans(dat[c("B-cells","Memory B-cells","naive B-cells","Class-switched memory B-cells","Plasma cells","pro B-cells"),])
  
  remove<-c("CD4+ memory T-cells","CD4+ naive T-cells",
            "CD4+ T-cells","CD4+ Tcm","CD4+ Tem",
            "CD8+ naive T-cells","CD8+ T-cells","CD8+ Tcm","CD8+ Tem",
            "aDC","cDC","DC","iDC","pDC",
            "Endothelial cells","ly Endothelial cells","mv Endothelial cells",
            "Macrophages","Macrophages M1","Macrophages M2",
            "B-cells","Memory B-cells","naive B-cells","Class-switched memory B-cells","pro B-cells",
            "Plasma cells","StromaScore","MicroenvironmentScore",
            "HSC","CLP","CMP","GMP","MEP","MPP",
            "Megakaryocytes","Erythrocytes","Platelets",
            "Smooth muscle","ImmuneScore")
  
  dat<-dat[which(!(rownames(dat)%in%remove)),]
  
  dat["CD4+ cells",]<-CD4
  dat["CD8+ cells",]<-CD8
  dat["Dendritic cells",]<-DC
  dat["Endothelial cells",]<-Endothelial
  dat["Macrophages",]<-Macrophage
  dat["B cells",]<-Bcell
  
  #Scale proportions so they sum to 1
  scaled_dat<-dat
  for (i in 1:ncol(dat)) {
    scaled_dat[,i]<-dat[,i]/sum(dat[,i])
  }
  return(scaled_dat)
}

celltypes<-fread("~/midway/xCell/CMC_GTEx_xCell_tpm_results.txt",data.table=FALSE)
rownames(celltypes)<-celltypes$V1
celltypes$V1<-NULL
celltypes<-prep_columns(celltypes)
celltypes<-as.data.frame(t(celltypes),stringsAsFactors=FALSE)
celltypes$prop_adipo<-celltypes$Adipocytes+celltypes$Preadipocytes
celltypes$sample<-rownames(celltypes)

#Get tissue IDs for each sample
tissues<-fread("~/midway/GTEx_Analysis_2017-06-05_v8_Annotations%2FGTEx_Analysis_2017-06-05_v8_Annotations_SampleAttributesDS.txt")
tissues<-tissues[,c("SAMPID","SMTSD")]
tissues$SMTSD<-gsub(" - ","_",tissues$SMTSD)
tissues$SMTSD<-gsub(" ","_",tissues$SMTSD)
tissues$SMTSD[tissues$SMTSD=="Cells_Cultured_fibroblasts"]<-"Cells_Transformed_fibroblasts"
colnames(tissues)<-c("sample","Tissue")
tissues<-rbind(tissues,cbind(colnames(celltypes)[!(grepl("GTEX",colnames(celltypes)))],rep("Brain_Cortex",length(cbind(colnames(celltypes)[!(grepl("GTEX",colnames(celltypes)))])))),use.names=FALSE)

#Merge tissue info with celltype estimations, keep only tissues also present in ldsc_seg
celltypes<-merge(celltypes,tissues,by="sample")
celltypes<-celltypes[celltypes$Tissue %in% dat$Tissue,]
celltypes<-celltypes %>% group_by(Tissue) %>% summarize(mean=mean(prop_adipo))
celltypes$tissue<-factor(celltypes$Tissue,levels=levels(dat$Tissue))
celltypes<-merge(celltypes,groups,by="Tissue")

#Plot mean by group
plot_dat<-celltypes %>% group_by(Group) %>% summarize(mean=mean(mean))
plot_dat$Group<-as.factor(plot_dat$Group)
A<-ggplot(plot_dat,aes(x=Group,y=mean))+
  geom_bar(aes(fill=Group),alpha=0.7,position="dodge",stat="identity")+
  theme_minimal()+
  scale_y_reverse("Proportion adipose cell types")+
  scale_fill_manual(values=rgb(pony_palette[,1:3]))+
  theme(axis.title=element_blank(),
        axis.text=element_text(size=20),
        legend.position = "none")+
  scale_x_discrete(breaks = NULL,limits=rev(levels(plot_dat$Group)))

A

pdf("~/medusa/papers/TWAS/intro/WHR_fat.pdf",width=5,height=3)
A
dev.off()

#Plot each tissue type
A<-ggplot(celltypes,aes(x=tissue,y=mean))+
  geom_bar(aes(fill=Group),alpha=0.7,position="dodge",stat="identity")+
  theme_minimal()+
  scale_y_continuous("Proportion adipose cell types")+
  scale_fill_manual(values=rgb(pony_palette[,1:3]))+
  theme(axis.title.y=element_text(size=15),
        axis.title.x=element_blank(),
        axis.text.x=element_text(size=10,angle=-90),
        legend.position="none")
A

pdf("~/medusa/papers/TWAS/supplement/intro/WHR_fat_all.pdf",width=7.5,height=4.5)
A
dev.off()
