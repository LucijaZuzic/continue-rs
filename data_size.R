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

oldw<-getOption('warn')
options(warn=-1)

# Data input: TEC and dTEC data derived from RINEX GPS observations taken at Darwin,NT,Australia
# using Gopi Seemala's GPS TEC software - days in 2014

iono1<-matrix(nrow=0,ncol=3)
geom1<-matrix(nrow=0,ncol=4)
dst1<-matrix(nrow=0,ncol=4)
all_previous_months<-0
timers<-c()
month_names<-c('January','February','March','April','May','June','July','August','September','October','November','December')
total_old<-0
total_new<-0
total_diff<-0
for(j in 1:12){ 
  month_file_name<-paste('darwin_data/data_all_',j,'.csv',sep='')
  data_month<-as.data.frame(read.csv(month_file_name,header=TRUE,sep=',',skip=0))
  old_len<-length(data_month$TEC)
  total_old<-total_old+old_len
  data_month<-data_month[data_month$TEC < limited,]
  new_len<-length(data_month$TEC)
  total_new<-total_new+new_len
  diff_len<-old_len-new_len
  total_diff<-total_diff+diff_len
  print(paste(month_names[j],old_len,new_len,diff_len,round(diff_len/old_len*100,2),round(new_len/old_len*100,2)))
  day_value<-all_previous_months + data_month['d']
  timer<-day_value + data_month['total.seconds']/86400
  timer<-round(timer,digits=4)
  diono<-as.data.frame(cbind(timer,data_month['TEC'],data_month['dTEC'],day_value,data_month['total.seconds'],data_month['m']))
  colnames(diono)<-c('time','TEC','dTEC','DOY','time_proper','month')
  iono1<-rbind(iono1,diono)
  geomi<-as.data.frame(cbind(timer,data_month['Bx'],data_month['By'],data_month['Bz']))
  colnames(geomi)<-c('time','Bx','By','Bz')
  geom1<-rbind(geom1,geomi)
  dsti<-as.data.frame(cbind(timer,data_month['Dst'],day_value))
  colnames(dsti)<-c('time','Dst','DOY')
  dst1<-rbind(dst1,dsti)
  timers<-rbind(timers,timer)
  all_previous_months<-all_previous_months + length(unique(data_month[['d']]))
}
print(paste(2014,total_old,total_new,total_diff,round(total_diff/total_old*100,2),round(total_new/total_old*100,2)))

# TEC & B data aggregation into single data frame per geomagnetic event

envi<-merge(iono1,geom1,by='time')

### Model development

Bx_prev<-c(envi$Bx[1], envi$Bx[-length(envi$Bx)])
By_prev<-c(envi$By[1], envi$By[-length(envi$By)])
Bz_prev<-c(envi$Bz[1], envi$Bz[-length(envi$Bz)])

envi$Bx_diff<-(envi$Bx - Bx_prev) / Bx_prev
envi$By_diff<-(envi$By - By_prev) / By_prev
envi$Bz_diff<-(envi$Bz - Bz_prev) / Bz_prev

envi<-cbind(envi['TEC'],envi$Bx,envi$By,envi$Bz,envi$Bx_diff,envi$By_diff,envi$Bz_diff)
colnames(envi)<-c('TEC','Bx','By','Bz','Bx_diff','By_diff','Bz_diff')
envi<-envi[-1,]

# Define a 80%/20% train/test split of the data

trainIndex<-createDataPartition(envi$TEC,p=0.80,list=FALSE)
dataTrain<-envi[trainIndex,]
dataTest<-envi[-trainIndex,]
print(paste('envi', length(envi$TEC)))
print(paste('dataTrain', length(dataTrain$TEC)))
print(paste('dataTest', length(dataTest$TEC)))

### END