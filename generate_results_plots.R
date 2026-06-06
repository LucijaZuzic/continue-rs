rm(list=ls())

# Declaration of the R libraries utilised

# Core data frameworks
library(data.table)
library(tidyverse)

# Miscellaneous UI
library(gridExtra)
library(import)

getCurrentFileLocation<-function()
{
  this_file<-commandArgs() %>% 
    tibble::enframe(name=NULL) %>%
    tidyr::separate(col=value,into=c('key','value'),sep='=',fill='right') %>%
    dplyr::filter(key=='--file') %>%
    dplyr::pull(value)
  if (length(this_file)==0)
  {      this_file<-rstudioapi::getSourceEditorContext()$path
  }
  return(dirname(this_file))
}

# Set path to the working directory containing:
# 1. INTERMAGNET data set-to be taken from the INTERMAGNET web-site,and reformatted 
# accordingly-# reformatted original data enclosed
# 2. TEC data set,as derived from the RINEX GPS observations using GPS TEC software-enclosed

setwd(getCurrentFileLocation())

limited<-150

# Synchronized arrays sorted logically from simplest baseline to heaviest ensemble

model_abbreviations<-c('LM','KNN','NN','RF','SGB')

dir.create(file.path(paste('TEC_',limited,sep='')),showWarnings=FALSE)
dir.create(file.path(paste('TEC_',limited,'/Klobuchar',sep='')),showWarnings=FALSE)

for (model_order in 1:length(model_abbreviations)) {
  dir.create(file.path(paste('TEC_',limited,'/',model_abbreviations[model_order],sep='')),showWarnings=FALSE)
}

oldw<-getOption('warn')
options(warn=-1)

# Data input: TEC and dTEC data derived from RINEX GPS observations taken at Darwin,NT,Australia
# using Gopi Seemala's GPS TEC software-days in 2014

### Model performance assessment

pdf(width=16,height=9,file=paste('TEC_',limited,'/Rotated_all_observed_predicted_',limited,'.pdf',sep=''))

par(mfrow=c(2,3),mar=c(5,5.5,4,1.5)+0.1)

for (model_order in 1:length(model_abbreviations)) {
  cart1_file<-as.data.frame(read.csv(paste('TEC_',limited,'/',model_abbreviations[model_order],'/',model_abbreviations[model_order],'_',limited,'.csv',sep=''),header=TRUE,sep=',',skip=0))
  plot(cart1_file$real,cart1_file$pred,cex.main=2,cex.lab=2,cex.axis=2,main=paste('Observed and ',model_abbreviations[model_order],' TEC [TECU]',sep=''),xlab='Observed TEC [TECU]',ylab=paste(model_abbreviations[model_order],'TEC [TECU]'))
  a<-c(0,50)
  b<-c(0,50)
  dabl<-as.data.frame(cbind(a,b))
  fit<-lm(b~a,data=dabl)
  abline(fit,col='red')
}

k_file<-as.data.frame(read.csv(paste('TEC_',limited,'/Klobuchar/Klobuchar_',limited,'.csv',sep=''),header=TRUE,sep=',',skip=0))
plot(k_file$real,k_file$pred,cex.main=2,cex.lab=2,cex.axis=2,main='Observed and Klobuchar model TEC [TECU]',xlab='Observed TEC [TECU]',ylab='Klobuchar model TEC [TECU]')
a<-c(0,50)
b<-c(0,50)
dabl<-as.data.frame(cbind(a,b))
fit<-lm(b~a,data=dabl)
abline(fit,col='red')

dev.off()