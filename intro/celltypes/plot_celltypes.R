#!/usr/bin/Rscript
library(data.table)
library(tidyverse)
library(yarrr) #For color funs
library(cowplot)
library(gridExtra)
setwd("~/papers/TWAS/TWAS_intro/")
pony_colors<-fread("~/medusa/papers/TWAS/pony_palette")

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

groups<-c("Adipose","Neural","Immune","Stromal","Immune"
          ,rep("Stromal",4),"Immune"
          ,rep("Stromal",2),"Immune","Stromal","Stromal",
          "Neural",rep("Immune",3),rep("Stromal",2),
          "Adipose","Stromal","Stromal",rep("Immune",7),
          "Stromal","Immune","Immune")

## Adipose
celltypes<-fread("~/midway/xCell/CMC_GTEx_xCell_tpm_results.txt",data.table=FALSE)
adipose_samples<-c("V1",scan("~/midway/expression/WHR/adipose_subcutaneous/adipose_subcutaneous_sample_IDs",what='character',sep='\n'))
celltypes_adipose<-celltypes[,colnames(celltypes)%in%adipose_samples]
rownames(celltypes_adipose)<-celltypes_adipose$V1
celltypes_adipose$V1<-NULL
celltypes_adipose<-prep_columns(celltypes_adipose)

adipose_dat<-as.data.frame(t(rbind(rownames(celltypes_adipose),
                                   rowMeans(celltypes_adipose),
                                   apply(celltypes_adipose,1,sd))),stringsAsFactors = FALSE)
colnames(adipose_dat)<-c("celltype","mean_prop","sd")
adipose_dat$mean_prop<-as.numeric(adipose_dat$mean_prop)
adipose_dat$sd<-as.numeric(adipose_dat$sd)
adipose_dat$groups<-groups
adipose_dat<-adipose_dat %>% arrange(groups,mean_prop)
adipose_dat$celltype<-factor(adipose_dat$celltype,levels=unique(adipose_dat$celltype))

##Cortex
celltypes_cortex<-celltypes[,!(grepl("GTEX",colnames(celltypes)))]
rownames(celltypes_cortex)<-celltypes_cortex$V1
celltypes_cortex$V1<-NULL
celltypes_cortex<-prep_columns(celltypes_cortex)

cortex_dat<-as.data.frame(t(rbind(rownames(celltypes_cortex),
                                  rowMeans(celltypes_cortex),
                                  apply(celltypes_cortex,1,sd))),stringsAsFactors = FALSE)
colnames(cortex_dat)<-c("celltype","mean_prop","sd")
cortex_dat$mean_prop<-as.numeric(cortex_dat$mean_prop)
cortex_dat$sd<-as.numeric(cortex_dat$sd)
cortex_dat$groups<-groups
cortex_dat<-cortex_dat[match(adipose_dat$celltype,cortex_dat$celltype),]
cortex_dat$celltype<-factor(cortex_dat$celltype,levels=unique(cortex_dat$celltype))

##Plot all celltypes
colors=c(rgb(pony_colors[5,1:3]),'gray70',rgb(pony_colors[8,1:3]),'gray70')

A<-ggplot(adipose_dat,aes(x=celltype,y=mean_prop))+
  geom_bar(aes(fill=groups),position="dodge",stat="identity")+
  geom_errorbar(aes(ymin=mean_prop-sd, ymax=mean_prop+sd), width=.2,
                position=position_dodge(.9))+
  theme_minimal()+
  scale_y_continuous("Estimated cell type proportion",limits=c(-0.02,0.2))+
  scale_fill_manual(values=c(rgb(pony_colors[5,1:3]),'gray70','gray70','gray70'))+
  theme(
    plot.title = element_text(size=24),
    axis.title=element_text(size=15),
    legend.position = "none",
    axis.title.x=element_blank(),
    panel.grid.major.x = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1,size=12))+
  ggtitle("Subcutaneous adipose RNA-seq")


C<-ggplot(cortex_dat,aes(x=celltype,y=mean_prop))+
  geom_bar(aes(fill=groups),position="dodge",stat="identity")+
  geom_errorbar(aes(ymin=mean_prop-sd, ymax=mean_prop+sd), width=.2,
                position=position_dodge(.9))+
  theme_minimal()+
  scale_y_continuous("Estimated cell type proportion",limits=c(-0.05,0.75))+
  scale_fill_manual(values=c('gray70','gray70',rgb(pony_colors[8,1:3]),'gray70'))+
  theme(
    plot.title = element_text(size=24),
    axis.title=element_text(size=15),
    legend.position = "none",
    axis.title.x=element_blank(),
    panel.grid.major.x = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1,size=12))+
  ggtitle("Frontal cortex RNA-seq")


pdf("~/medusa/papers/TWAS/supplement/intro/celltypes_all.pdf",width=7,height=7)
grid.arrange(A,C,nrow=2)
dev.off()


######################### Plot by sex
# Get sex information
#celltypes_tissue<-celltypes_tissue[,grepl("GTEX",colnames(celltypes_tissue))]
sex<-fread("~/midway/GTEx_Analysis_2017-06-05_v8_Annotations%2FGTEx_Analysis_2017-06-05_v8_Annotations_SubjectPhenotypesDS.txt")
sex<-sex[,c('SUBJID',"SEX")]
subjs<-paste("GTEX",sapply(strsplit(colnames(celltypes_adipose),'-'),'[[',2),sep='-')
sex<-sex$SEX[match(subjs,sex$SUBJID)] #In same order as GTEx samples

props_F<-as.data.frame(t(rbind(rownames(celltypes_adipose),
                               rowMeans(celltypes_adipose[,sex==2]),
                               apply(celltypes_adipose[,sex==2],1,sd))),stringsAsFactors = FALSE)
colnames(props_F)<-c("celltype","mean_prop","sd")
props_F$mean_prop<-as.numeric(props_F$mean_prop)
props_F$sd<-as.numeric(props_F$sd)
props_F$celltype<-factor(props_F$celltype,levels=unique(props_F$celltype))
props_F$sex="F"

props_M<-as.data.frame(t(rbind(rownames(celltypes_adipose),
                               rowMeans(celltypes_adipose[,sex==1]),
                               apply(celltypes_adipose[,sex==1],1,sd))),stringsAsFactors = FALSE)
colnames(props_M)<-c("celltype","mean_prop","sd")
props_M$mean_prop<-as.numeric(props_M$mean_prop)
props_M$sd<-as.numeric(props_M$sd)
props_M$celltype<-factor(props_M$celltype,levels=unique(props_M$celltype))
props_M$sex="M"

props<-rbind(props_F,props_M)
props$celltype<-factor(props$celltype,levels=adipose_dat$celltype)
pdf("~/midway/xCell/adipose_subcutaneous_xCell_by_sex.pdf",width=7,height=4)
C<-ggplot(props,aes(x=celltype,y=mean_prop,fill=sex))+
  geom_bar(position="dodge",stat="identity")+
  geom_errorbar(aes(ymin=mean_prop-sd, ymax=mean_prop+sd), width=.2,
                position=position_dodge(0.9))+
  theme_minimal()+
  scale_y_continuous("Estimated cell type proportion",limits=c(-0.05,0.75))+
  scale_fill_manual(values=c(rgb(pony_colors[14,1:3]),rgb(pony_colors[10,1:3])))+
  theme(
    plot.title = element_text(size=24),
    axis.title=element_text(size=15),
    legend.position = c(0.75, 0.68),
    legend.title=element_blank(),
    axis.title.x=element_blank(),
    panel.grid.major.x = element_blank(),
    axis.text.y=element_text(size=12),
    axis.text.x = element_text(angle = 45, hjust = 1,size=12))+
  ggtitle(paste("Subcutaneous adipose estimated \n cell type proportions by sex",sep=''))
C
dev.off()
