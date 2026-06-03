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

kl.pr<-c(0.0864608287811279,14.4241651168816,0.446739923128779,11.6170828925552)
lm.pr<-c(0.0471348762512207,14.0015329,0.3612796,11.0134747)
cart1.pr<-c(2.56556634902954*60,12.1051081,0.5225813,8.8273016)
gamb.pr<-c(9.31457554896673*60,11.9172233,0.5378511,8.8369546)
bcart.pr<-c(18.0784080306689*60,12.0472057,0.5271387,8.7850788)
sgb.pr<-c(21.2666043361028*60,11.0682894,0.6013171,7.9176253)

### Model performance assessment

pdf(file=paste('TEC_',limited,'/MAE_RMSE_R2_Time_',limited,'.pdf',sep=''))
par(mfrow=c(2,2))

# Create the data for the MAE chart

H<-c(lm.pr[4], cart1.pr[4],sgb.pr[4],bcart.pr[4],gamb.pr[4],kl.pr[4])
M<-c('LM','CART1','SGB','BCART','GAMB','K')

# Plot the bar chart

bar_centers<-barplot(H,cex.names=0.7,names.arg=M,xlab='Method',ylab='MAE [TECU]',ylim=c(0,max(H)*1.15),col='gray',main='MAE [TECU]',border='red')
text(x=bar_centers,y=H,labels=round(H,2),pos=3,cex=0.9,col="black")

# Create the data for the RMSE chart

H<-c(lm.pr[2],cart1.pr[2],sgb.pr[2],bcart.pr[2],gamb.pr[2],kl.pr[2])
M<-c('LM','CART1','SGB','BCART','GAMB','K')

# Plot the bar chart

bar_centers<-barplot(H,cex.names=0.7,names.arg=M,xlab='Method',ylab='RMSE [TECU]',ylim=c(0,max(H)*1.15),col='gray',main='RMSE [TECU]',border='red')
text(x=bar_centers,y=H,labels=round(H,2),pos=3,cex=0.9,col="black")

# Create the data for the adjR2 chart

H<-c(lm.pr[3]*100,cart1.pr[3]*100,sgb.pr[3]*100,bcart.pr[3]*100,gamb.pr[3]*100,kl.pr[3]*100)
M<-c('LM','CART1','SGB','BCART','GAMB','K')

# Plot the bar chart

bar_centers<-barplot(H,cex.names=0.7,names.arg=M,xlab='Method',ylab='Adjusted R2 [%]',ylim=c(0,max(H)*1.15),col='grey',main='Adjusted R2 [%]',border='red')
text(x=bar_centers,y=H,labels=round(H,2),pos=3,cex=0.9,col="black")

# Create the data for the model development time chart

H<-c(lm.pr[1],cart1.pr[1],sgb.pr[1],bcart.pr[1],gamb.pr[1],kl.pr[1])
M<-c('LM','CART1','SGB','BCART','GAMB','K')

# Plot the bar chart

bar_centers<-barplot(H,cex.names=0.7,names.arg=M,xlab='Method',ylab='Time elapsed [s]',ylim=c(0,max(H)*1.15),col='gray',main='Model development time [s]',border='red')
text(x=bar_centers,y=H,labels=round(H,2),pos=3,cex=0.9,col="black")
dev.off()