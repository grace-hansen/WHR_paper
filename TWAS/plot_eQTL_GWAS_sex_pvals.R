#!/usr/bin/Rscript
library(data.table)
library(tidyverse)
library(eulerr)
library(yarrr) #For color funs
library(optparse)
library(gridExtra)
#arguments <- parse_args(Optionparser(usage = "plot_eQTL_GWAS_sex.R ",option_list=list()),
#                        positional_arguments = 2)
#opt=arguments$opt
trait<-"WHR"
tissue<-"adipose_subcutaneous"
data_source<-"GTEx_v8"
sex<-"F"
pony_colors<-fread("~/medusa/papers/TWAS/pony_palette")
setwd("~/medusa/papers/TWAS/TWAS/")

############################ For color manipulation ############################
darken <- function(color, factor=1.2){
  col <- col2rgb(color)
  col <- col/factor
  col <- rgb(t(col), maxColorValue=255)
  col
}
################################################################################

#Grab significant genes
genes<-fread(paste("~/midway/expression/",trait,"_",sex,"/",tissue,"/results/posthoc/sig_genes",sep=''))

#Grab top eQTL variant for each gene
TWAS<-fread(paste("~/midway/expression/",trait,"_",sex,"/",tissue,"/results/GTEx_v8.all.dat.top",sep=''))

#Grab p-value values of eQTL for each gene in each sex
get_vals<-function() {
  eQTL_p_F<-numeric(nrow(TWAS))
  eQTL_p_M<-numeric(nrow(TWAS))
  for (i in 1:nrow(TWAS)) {
    chr<-TWAS$CHR[i]
    gene<-TWAS$ID[i]
    var<-TWAS$EQTL.ID[i]
    p_F<-system(paste("awk -F ' ' '$1==\"",gene,"\"&&$2==\"",var,"\" {print $0}' ~/midway/QTL_analyses/eQTL/GTEx_F/results/chr",chr,"_eQTL.txt | cut -f4 -d' '",sep=''),intern=TRUE)
    p_M<-system(paste("awk -F ' ' '$1==\"",gene,"\"&&$2==\"",var,"\" {print $0}' ~/midway/QTL_analyses/eQTL/GTEx_M/results/chr",chr,"_eQTL.txt | cut -f4 -d' '",sep=''),intern=TRUE)
    if (length(p_F)>0 && length(p_M)>0) {
      eQTL_p_F[i]<-p_F
      eQTL_p_M[i]<-p_M
    }
  }
  #Grab p-values of eQTL variant for GWAS in each sex
  GWAS_p_F<-numeric(nrow(TWAS))
  GWAS_p_M<-numeric(nrow(TWAS))
  for (i in 1:nrow(TWAS)) {
    chr<-TWAS$CHR[i]
    gene<-TWAS$ID[i]
    var<-TWAS$EQTL.ID[i]
    p_F<-system(paste("zcat ~/midway/expression/WHR_F/Whradjbmi.giant-ukbb.meta-analysis.females.23May2018.HapMap2_only.txt.gz | grep ",var,": | cut -f9 -d' '",sep=''),intern=TRUE)
    p_M<-system(paste("zcat ~/midway/expression/WHR_M/Whradjbmi.giant-ukbb.meta-analysis.males.23May2018.HapMap2_only.txt.gz | grep ",var,": | cut -f9 -d' '",sep=''),intern=TRUE)
    if (length(p_F)>0 && length(p_M)>0) {
      GWAS_p_F[i]<-p_F
      GWAS_p_M[i]<-p_M
    }
  }
  return(list(eQTL_p_F,eQTL_p_M,GWAS_p_F,GWAS_p_M))
}
#vals<-get_vals()
#eQTL_p_F<-vals[[1]]
#eQTL_p_M<-vals[[2]]
#GWAS_p_F<-vals[[3]]
#GWAS_p_M<-vals[[4]]

#Load p-values so I don't have to find them again
plot_dat<-fread("eQTL_GWAS_pvals.txt")

#plot
#plot_dat<-cbind(TWAS[,c("ID","EQTL.ID")],eQTL_p_F,eQTL_p_M,GWAS_p_F,GWAS_p_M)
#write.table(plot_dat,"eQTL_GWAS_pvals.txt",quote=FALSE,row.names=FALSE,sep='\t')
plot_dat<-as.data.frame(plot_dat[rowSums(plot_dat!=0)==6,])
for (i in 3:6) {
  plot_dat[,i]<-as.numeric(plot_dat[,i])
}

eF<--log10(plot_dat$eQTL_p_F)
eM<--log10(plot_dat$eQTL_p_M)
slope_E<-lm(eM ~ eF)$coefficients[2]
int_E<-lm(eM ~ eF)$coefficients[1]
E<-ggplot(plot_dat,aes(x=-log10(eQTL_p_F),y=-log10(eQTL_p_M)))+
  geom_abline(linetype="dotted")+
  geom_abline(slope=slope_E,intercept=int_E,size=1)+
  geom_point(size=2,color=darken(rgb(pony_colors[3,1:3])))+
  annotate("text",x=12,y=25,label=paste("slope=",round(slope_E,digits=4),sep=''),size=5)+
  scale_x_continuous("-log10 eQTL p-value, females",limits=c(0,40))+
  scale_y_continuous("-log10 eQTL p-value, males",limits=c(0,40))+
  theme_minimal()+
  theme(axis.text = element_text(size=16),
        axis.title = element_text(size=16))

gF<--log10(plot_dat$GWAS_p_F)
gM<--log10(plot_dat$GWAS_p_M)
slope_G<-lm(gM ~ gF)$coefficients[2]
int_G<-lm(gM ~ gF)$coefficients[1]
G<-ggplot(plot_dat,aes(x=-log10(GWAS_p_F),y=-log10(GWAS_p_M)))+
  geom_abline(linetype="dotted")+
  geom_abline(slope=slope_G,intercept=int_G,size=1)+
  geom_point(size=2,color=rgb(pony_colors[7,1:3]))+
  annotate("text",x=14,y=37,label=paste("slope=",round(slope_G,digits=4),sep=''),size=5)+
  scale_x_continuous("-log10 GWAS p-value, females",limits=c(0,60))+
  scale_y_continuous("-log10 GWAS p-value, males",limits=c(0,60))+
  theme_minimal()+
  theme(axis.text = element_text(size=16),
        axis.title = element_text(size=16))
G
pdf("eQTL_GWAS_sex_pvals_scatter.pdf",width=4,height=8)
grid.arrange(E,G,nrow=2)
dev.off()

eQTL<-pivot_longer(plot_dat,c(eQTL_p_F,eQTL_p_M),names_to="sex",values_to="eQTL")
GWAS<-pivot_longer(plot_dat,c(GWAS_p_F,GWAS_p_M),names_to="sex",values_to="GWAS")
eQTL$sex[eQTL$sex=="eQTL_p_F"]<-"F"
eQTL$sex[eQTL$sex=="eQTL_p_M"]<-"M"
GWAS$sex[GWAS$sex=="GWAS_p_F"]<-"F"
GWAS$sex[GWAS$sex=="GWAS_p_M"]<-"M"
plot_dat<-eQTL
plot_dat$GWAS<-GWAS$GWAS
plot_dat$GWAS_p_F<-NULL
plot_dat$GWAS_p_M<-NULL

E<-ggplot(plot_dat,aes(x=sex,y=-log10(eQTL),group=EQTL.ID))+
  geom_violin(aes(x=sex,y=-log10(eQTL),group=sex),fill=rgb(pony_colors[3,1:3]),draw_quantiles=c(0.5))+
  geom_point()+
  geom_line(alpha=0.2)+
  scale_y_continuous("-log10 eQTL p-value",limits=c(0,50))+
  theme_minimal()+
  theme(axis.text = element_text(size=20),
        axis.title = element_text(size=20))
E

G<-ggplot(plot_dat,aes(x=sex,y=-log10(GWAS),group=EQTL.ID))+
  geom_violin(aes(x=sex,y=-log10(GWAS),group=sex),fill=rgb(pony_colors[7,1:3]),draw_quantiles=c(0.5))+
  geom_point()+
  geom_line(alpha=0.2)+
  scale_y_continuous("-log10 GWAS p-value",limits=c(0,70))+
  theme_minimal()+
  theme(axis.text = element_text(size=20),
        axis.title = element_text(size=20))

grid.arrange(E,G,nrow=1)

pdf("eQTL_GWAS_sex_pvals_violin.pdf",width=8,height=6)
grid.arrange(E,G,nrow=1)
dev.off()
