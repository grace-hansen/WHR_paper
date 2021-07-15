#!/usr/bin/Rscript
library(data.table)
library(tidyverse)
library(optparse)
library(gridExtra)
pony_colors<-fread("~/medusa/papers/TWAS/pony_palette")
arguments <- parse_args(OptionParser(usage = "plot_gen_corr.R <trait1> <trait2>",option_list=list()),
                        positional_arguments = 2)
opt=arguments$opt
traits<-c(arguments$args[1],arguments$args[2])
##################### Author: Grace Hansen ###################
# The purpose of this script is to plot h2 and genetic correlation estimates for each sex and between each sex from ldsc. 

##################### Genetic correlation ####################

dat<-as.data.frame(matrix(ncol=3,nrow=2),stringsAsFactors=FALSE)
colnames(dat)<-c("trait","corr","se")

for (i in 1:length(traits)) {
  dat$trait[i]<-traits[i]
  dat$corr[i]<-system(paste("grep -m 2 Correlation ~/midway/ldsc/",traits[i],"_F_",traits[i],"_M_corr/",traits[i],"_F_",traits[i],"_M.log | tail -1 | cut -f3 -d' '",sep=''),intern=TRUE)
  dat$se[i]<-system(paste("grep -m 2 Correlation ~/midway/ldsc/",traits[i],"_F_",traits[i],"_M_corr/",traits[i],"_F_",traits[i],"_M.log | tail -1 | cut -f4 -d' ' | cut -f2 -d'(' | cut -f1 -d')'",sep=''),intern=TRUE)
}
dat$corr<-as.numeric(dat$corr)
dat$se<-as.numeric(dat$se)

##Plot
dat$trait[which(dat$trait=="WHR")]<-"WHRadjBMI" #Rename WHR to WHRadjBMI if necessary
C<-ggplot(dat)+
  geom_col(aes(x=trait,y=corr,fill=trait))+
  geom_errorbar(aes(x=trait,ymin=corr-se,ymax=corr+se),width=0.4,size=1.25)+
  ggtitle("Genetic correlation between \nmale and female genetic risk")+
  coord_cartesian(ylim=c(0,1))+
  theme_minimal()+
  scale_y_continuous(name="Correlation (r2)")+
  scale_fill_manual(values=c(rgb(pony_colors[3,1:3]),rgb(pony_colors[10,1:3])))+
  theme(plot.title=element_blank(),
        legend.position="none",
        axis.text.y=element_blank(),
        axis.text.x=element_text(size=40),
        axis.title.x=element_text(size=44),
        axis.title.y=element_blank())+
  coord_flip()

pdf(paste("~/midway/ldsc/gen_corr_sex_",traits[1],"_",traits[2],".pdf",sep=''),width=7,height=4)
C
dev.off()

################## Heritability, specific #####################

for (i in 1:length(traits)) {
  dat<-as.data.frame(matrix(ncol=4,nrow=2),stringsAsFactors=FALSE)
  colnames(dat)<-c("trait","sex","h2","se")
  dat$trait[1]<-traits[i]
  dat$sex[1]<-"F"
  dat$h2[1]<-system(paste("grep Observed ~/midway/ldsc/",traits[i],"_F/",traits[i],"_F_h2.log | cut -f5 -d' '",sep=''),intern=TRUE)
  dat$se[1]<-system(paste("grep Observed ~/midway/ldsc/",traits[i],"_F/",traits[i],"_F_h2.log | cut -f6 -d' ' | cut -f2 -d'(' | cut -f1 -d')'",sep=''),intern=TRUE)
  dat$trait[2]<-traits[i]
  dat$sex[2]<-"M"
  dat$h2[2]<-system(paste("grep Observed ~/midway/ldsc/",traits[i],"_M/",traits[i],"_M_h2.log | cut -f5 -d' '",sep=''),intern=TRUE)
  dat$se[2]<-system(paste("grep Observed ~/midway/ldsc/",traits[i],"_M/",traits[i],"_M_h2.log | cut -f6 -d' ' | cut -f2 -d'(' | cut -f1 -d')'",sep=''),intern=TRUE)
  
  dat$h2<-as.numeric(dat$h2)
  dat$se<-as.numeric(dat$se)
  
  if (traits[i]=="obesity") {
    col=rgb(pony_colors[3,1:3])
  } else if (traits[i]=="WHR") {
    col=rgb(pony_colors[10,1:3])
  }
  
  ##Plot
  C<-ggplot(dat)+
    geom_col(aes(x=sex,y=h2,fill=trait))+
    geom_errorbar(aes(x=sex,ymin=h2-se,ymax=h2+se),width=0.4,size=1.25)+
    scale_y_continuous(name="Heritability (h2)")+
    scale_x_discrete(name="Sex")+
    theme_minimal()+
    ggtitle(paste(traits[i]," heritability by sex",sep=''))+
    scale_fill_manual(values=c(col))+
    theme(plot.title=element_blank(),
          legend.position="none",
          axis.text=element_text(size=40),
          axis.title.x=element_text(size=44),
          axis.title.y=element_blank())
  C
  
  pdf(paste("~/midway/ldsc/h2_sex_",traits[i],".pdf",sep=''),width=6,height=10)
  C
  dev.off()
}

######### Heritability, sex-combined ##################

dat<-as.data.frame(matrix(ncol=3,nrow=2),stringsAsFactors=FALSE)
colnames(dat)<-c("trait","h2","se")
dat$trait[1]<-traits[1]
dat$h2[1]<-system(paste("grep Observed ~/midway/ldsc/",traits[1],"/",traits[1],"_h2.log | cut -f5 -d' '",sep=''),intern=TRUE)
dat$se[1]<-system(paste("grep Observed ~/midway/ldsc/",traits[1],"/",traits[1],"_h2.log | cut -f6 -d' ' | cut -f2 -d'(' | cut -f1 -d')'",sep=''),intern=TRUE)
dat$trait[2]<-traits[2]
dat$h2[2]<-system(paste("grep Observed ~/midway/ldsc/",traits[2],"/",traits[2],"_h2.log | cut -f5 -d' '",sep=''),intern=TRUE)
dat$se[2]<-system(paste("grep Observed ~/midway/ldsc/",traits[2],"/",traits[2],"_h2.log | cut -f6 -d' ' | cut -f2 -d'(' | cut -f1 -d')'",sep=''),intern=TRUE)
  
dat$h2<-as.numeric(dat$h2)
dat$se<-as.numeric(dat$se)
  
col<-c(rgb(pony_colors[3,1:3]),rgb(pony_colors[10,1:3]))
       
  ##Plot
C<-ggplot(dat)+
  geom_col(aes(x=trait,y=h2,fill=trait))+
  geom_errorbar(aes(x=trait,ymin=h2-se,ymax=h2+se),width=0.4,size=1.25)+
  scale_y_continuous(name="Heritability (h2)")+
  scale_x_discrete(name="Sex")+
  theme_minimal()+
  scale_fill_manual(values=c(col))+
  theme(plot.title=element_blank(),
        legend.position="none",
        axis.text=element_text(size=40),
        axis.title=element_text(size=44))
C
  
pdf(paste("~/midway/ldsc/h2_",paste(traits,collapse = "_"),".pdf",sep=''),width=6,height=6)
C
dev.off()
  