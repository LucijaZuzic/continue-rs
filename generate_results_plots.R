rm(list=ls())

# Declaration of the R libraries utilised

library(data.table)
library(AppliedPredictiveModeling)
library(DataExplorer)
library(corrplot)
library(glmnet)
library(caret)
library(car)
library(forecast)
library(relaimpo)
library(DALEX)
library(monomvn)
library(mboost)
library(import)
library(plyr)
library(monomvn)
library(gridExtra)
library(tidyverse)
library(fitdistrplus)

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
# 1. INTERMAGNET data set - to be taken from the INTERMAGNET web-site,and reformatted 
# accordingly - # reformatted original data enclosed
# 2. TEC data set,as derived from the RINEX GPS observations using GPS TEC software - enclosed

setwd(getCurrentFileLocation())

limited<-150

dir.create(file.path(paste('TEC_',limited,sep='')),showWarnings=FALSE)
dir.create(file.path(paste('TEC_',limited,'/Klobuchar',sep='')),showWarnings=FALSE)
dir.create(file.path(paste('TEC_',limited,'/GAMB',sep='')),showWarnings=FALSE)
dir.create(file.path(paste('TEC_',limited,'/BCART',sep='')),showWarnings=FALSE)
dir.create(file.path(paste('TEC_',limited,'/CART1',sep='')),showWarnings=FALSE)
dir.create(file.path(paste('TEC_',limited,'/LM',sep='')),showWarnings=FALSE)
dir.create(file.path(paste('TEC_',limited,'/SGB',sep='')),showWarnings=FALSE)

### Model performance assessment

pdf(width=16,height=9,file=paste('TEC_',limited,'/Rotated_all_observed_predicted_',limited,'.pdf',sep=''))
par(mfrow=c(2,3),mar=c(5,5.5,4,1.5)+0.1)
lm_file<-as.data.frame(read.csv(paste('TEC_',limited,'/LM/LM_',limited,'.csv',sep=''),header=TRUE,sep=',',skip=0))
plot(lm_file$real,lm_file$pred,cex.main=2,cex.lab=2,cex.axis=2,main='Observed and linear regression TEC [TECU]',xlab='Observed TEC [TECU]',ylab='Linear regression TEC [TECU]')
a<-c(0,50)
b<-c(0,50)
dabl<-as.data.frame(cbind(a,b))
fit<-lm(b~a,data=dabl)
abline(fit,col='red')

cart1_file<-as.data.frame(read.csv(paste('TEC_',limited,'/CART1/CART1_',limited,'.csv',sep=''),header=TRUE,sep=',',skip=0))
plot(cart1_file$real,cart1_file$pred,cex.main=2,cex.lab=2,cex.axis=2,main='Observed and CART1 TEC [TECU]',xlab='Observed TEC [TECU]',ylab='CART1 TEC [TECU]')
a<-c(0,50)
b<-c(0,50)
dabl<-as.data.frame(cbind(a,b))
fit<-lm(b~a,data=dabl)
abline(fit,col='red')

bcart_file<-as.data.frame(read.csv(paste('TEC_',limited,'/BCART/BCART_',limited,'.csv',sep=''),header=TRUE,sep=',',skip=0))
plot(bcart_file$real,bcart_file$pred,cex.main=2,cex.lab=2,cex.axis=2,main='Observed and BCART TEC [TECU]',xlab='Observed TEC [TECU]',ylab='BCART TEC [TECU]')
a<-c(0,50)
b<-c(0,50)
dabl<-as.data.frame(cbind(a,b))
fit<-lm(b~a,data=dabl)
abline(fit,col='red')

gamb_file<-as.data.frame(read.csv(paste('TEC_',limited,'/GAMB/GAMB_',limited,'.csv',sep=''),header=TRUE,sep=',',skip=0))
plot(gamb_file$real,gamb_file$pred,cex.main=2,cex.lab=2,cex.axis=2,main='Observed and GAMB TEC [TECU]',xlab='Observed TEC [TECU]',ylab='GAMB TEC [TECU]')
a<-c(0,50)
b<-c(0,50)
dabl<-as.data.frame(cbind(a,b))
fit<-lm(b~a,data=dabl)
abline(fit,col='red')

sgb_file<-as.data.frame(read.csv(paste('TEC_',limited,'/SGB/SGB_',limited,'.csv',sep=''),header=TRUE,sep=',',skip=0))
plot(sgb_file$real,sgb_file$pred,cex.main=2,cex.lab=2,cex.axis=2,main='Observed and SGB TEC [TECU]',xlab='Observed TEC [TECU]',ylab='SGB TEC [TECU]')
a<-c(0,50)
b<-c(0,50)
dabl<-as.data.frame(cbind(a,b))
fit<-lm(b~a,data=dabl)
abline(fit,col='red')

k_file<-as.data.frame(read.csv(paste('TEC_',limited,'/Klobuchar/Klobuchar_',limited,'.csv',sep=''),header=TRUE,sep=',',skip=0))
plot(k_file$real,k_file$pred,cex.main=2,cex.lab=2,cex.axis=2,main='Observed and Klobuchar model TEC [TECU]',xlab='Observed TEC [TECU]',ylab='Klobuchar model TEC [TECU]')
a<-c(0,50)
b<-c(0,50)
dabl<-as.data.frame(cbind(a,b))
fit<-lm(b~a,data=dabl)
abline(fit,col='red')
dev.off()