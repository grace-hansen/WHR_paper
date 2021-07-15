#!/usr/bin/Rscript
library(data.table)
library(tidyverse)
library(pracma)
library(ggrepel)


setwd("~/medusa/papers/TWAS/TWAS/locus_zoom/")
sex<-args[1]
chr<-7
rsid<-"rs1534696"
setwd(sex)

######### Author: Grace Hansen #########
#This script makes sex-specific locus zoom plots near a given variant from GTEx ld data and sex-specific and sex-combined WHRadjBMI GWAS data. 

#Collect data
if (sex=='F') {
  cmd=paste("plink --bfile /home/grace/midway/LD/GTEx_v8_females/",chr," --ld-snp ",rsid," --ld-window 200000 --ld-window-r2 0 --maf 0.01 --out F_LD --r2",sep='')
  system(cmd)  
} else if (sex == "M") {
  cmd=paste("plink --bfile /home/grace/midway/LD/GTEx_v8_males/",chr," --ld-snp ",rsid," --ld-window 200000 --ld-window-r2 0 --maf 0.01 --out M_LD --r2",sep='')
  system(cmd)  
}

#Load data
ld<-fread(paste(sex,"_LD.ld",sep=''))
snp_pos<-ld$BP_A[1]
ld$rsid<-ld$SNP_B
ld$pos<-ld$BP_B
ld$chr<-ld$CHR_B
ld<-ld[,c("chr","pos","rsid","R2")]
if (sex=="F") {
  WHR_sex=fread(cmd="zcat ~/midway/expression/WHR_F/Whradjbmi.giant-ukbb.meta-analysis.females.23May2018.HapMap2_only_hg38.txt.gz")
} else if (sex=="M") {
  WHR_sex=fread(cmd="zcat  ~/midway/expression/WHR_M/Whradjbmi.giant-ukbb.meta-analysis.males.23May2018.HapMap2_only_hg38.txt.gz")
}
WHR=fread(cmd="zcat ~/midway/expression/WHR/Whradjbmi.giant-ukbb.meta-analysis.combined.23May2018.HapMap2_only_hg38.txt.gz")
WHR_sex$rsid<-sapply(strsplit(WHR_sex$SNP,':'),'[[',1)
WHR$rsid<-sapply(strsplit(WHR$SNP,':'),'[[',1)

dat<-merge(ld,WHR,by="rsid")
dat<-merge(dat,WHR_sex,by="rsid",suffixes=c("",paste("_",sex,sep='')))
dat$LD=character()
dat$LD[dat$R2<=0.2]="0-0.2"
dat$LD[dat$R2>0.2 & dat$R2<=0.4]="0.2-0.4"
dat$LD[dat$R2>0.4 & dat$R2<=0.6]="0.4-0.6"
dat$LD[dat$R2>0.6 & dat$R2<=0.8]="0.6-0.8"
dat$LD[dat$R2>0.8]="0.8-1.0"

dat<-dat[abs(dat$pos-snp_pos)<=200000,]
dat$pos=dat$pos/1000000

#Plot data
dat<-as.data.frame(dat)
var<-rsid
if (sex=='F') {
  G<-ggplot(dat)+
    geom_area(aes(x=pos,y=-log10(P),fill="Sex-combined p value"))+
    geom_point(aes(x=pos,y=-log10(P_F),color=LD),size=2)+
    geom_point(aes(x=dat$pos[dat$rsid==var],y=-log10(dat$P_F[dat$rsid==var])),color="firebrick3",size=2)+
    scale_y_continuous(limits=c(0,70))+
    scale_x_continuous(name="Position (MB)")+
    scale_color_manual(limits=c("0-0.2","0.2-0.4","0.4-0.6","0.6-0.8","0.8-1.0"),values = c("dodgerblue4","dodgerblue1","green2","goldenrod2","firebrick3"))+
    scale_fill_manual(values="gray80")+
    guides(color = guide_legend(override.aes = list(size=5)))+
    theme_minimal()+
    theme(panel.border = element_blank(), axis.line.x = element_line(),
          panel.grid.major.x = element_blank(),
          panel.grid.minor.x = element_blank(),
          legend.title=element_text(size=16),
          legend.text=element_text(size=16),
          axis.text=element_text(size=22),
          axis.title=element_text(size=22))
  G
} else if (sex == "M") {
  G<-ggplot(dat)+
    geom_area(aes(x=pos,y=-log10(P),fill="Sex-combined p value"))+
    geom_point(aes(x=pos,y=-log10(P_M),color=LD),size=2)+
    geom_point(aes(x=dat$pos[dat$rsid==var],y=-log10(dat$P_M[dat$rsid==var])),color="firebrick3",size=2)+
    scale_y_continuous(limits=c(0,70))+
    scale_x_continuous(name="Position (MB)")+
    scale_color_manual(limits=c("0-0.2","0.2-0.4","0.4-0.6","0.6-0.8","0.8-1.0"),values = c("dodgerblue4","dodgerblue1","green2","goldenrod2","firebrick3"))+
    scale_fill_manual(values="gray80")+
    guides(color = guide_legend(override.aes = list(size=5)))+
    theme_minimal()+
    theme(panel.border = element_blank(), axis.line.x = element_line(),
          panel.grid.major.x = element_blank(),
          panel.grid.minor.x = element_blank(),
          legend.title=element_text(size=16),
          legend.text=element_text(size=16),
          axis.text=element_text(size=22),
          axis.title=element_text(size=22))
  G
}
pdf(paste("~/medusa/papers/TWAS/TWAS/locus_zoom/",rsid,"_",sex,"_locuszoom.pdf",sep=''),width=10,height=5)
G
dev.off()
