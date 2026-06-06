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

dir.create(file.path(paste('TEC_',limited,sep='')),showWarnings=FALSE)

oldw<-getOption('warn')
options(warn=-1)

# Data input: TEC and dTEC data derived from RINEX GPS observations taken at Darwin,NT,Australia
# using Gopi Seemala's GPS TEC software-days in 2014

final_metrics_df<-as.data.frame(read.csv(paste('TEC_',limited,'/All_Models_Metrics_',limited,'.csv',sep=''),header=TRUE,sep=',',skip=0))
final_metrics_df<-final_metrics_df %>% replace(is.na(.), 0)

### Model performance assessment

pdf(file=paste('TEC_',limited,'/MAE_RMSE_R2_Time_',limited,'.pdf',sep=''))
par(mfrow=c(2,2))

# Create the data for the MAE chart

H<-as.numeric(final_metrics_df$MAE)
M<-final_metrics_df$name

# Plot the bar chart

bar_centers<-barplot(H,cex.names=0.9,las=2,names.arg=M,xlab='Method',ylab='MAE [TECU]',ylim=c(0,max(H)*1.4),col='gray',main='MAE [TECU]',border='red')
text(x=bar_centers,y=H,labels=round(H,2),adj=c(-0.25,0.5),srt=90,cex=0.9,col="black")

# Create the data for the RMSE chart

H<-as.numeric(final_metrics_df$RMSE)
M<-final_metrics_df$name

# Plot the bar chart

bar_centers<-barplot(H,cex.names=0.9,las=2,names.arg=M,xlab='Method',ylab='RMSE [TECU]',ylim=c(0,max(H)*1.4),col='gray',main='RMSE [TECU]',border='red')
text(x=bar_centers,y=H,labels=round(H,2),adj=c(-0.25,0.5),srt=90,cex=0.9,col="black")

# Create the data for the adjR2 chart

H<-as.numeric(final_metrics_df$adjRsquared)*100
M<-final_metrics_df$name

# Plot the bar chart

bar_centers<-barplot(H,cex.names=0.9,las=2,names.arg=M,xlab='Method',ylab='Adjusted R2 [%]',ylim=c(0,max(H)*1.4),col='grey',main='Adjusted R2 [%]',border='red')
text(x=bar_centers,y=H,labels=round(H,2),adj=c(-0.25,0.5),srt=90,cex=0.9,col="black")

# Create the data for the model development time chart

H<-as.numeric(final_metrics_df$training)
M<-final_metrics_df$name

# Plot the bar chart

bar_centers<-barplot(H,cex.names=0.9,las=2,names.arg=M,xlab='Method',ylab='Model development time [s]',ylim=c(0,max(H)*1.4),col='gray',main='Model development time [s]',border='red')
text(x=bar_centers,y=H,labels=round(H,2),adj=c(-0.25,0.5),srt=90,cex=0.9,col="black")
dev.off()