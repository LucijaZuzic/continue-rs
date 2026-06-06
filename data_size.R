rm(list=ls())

# Declaration of the R libraries utilised

# Core data frameworks
library(data.table)
library(tidyverse)

# Model wrapper
library(caret)

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

oldw<-getOption('warn')
options(warn=-1)

# Data input: TEC and dTEC data derived from RINEX GPS observations taken at Darwin,NT,Australia
# using Gopi Seemala's GPS TEC software-days in 2014

iono1<-matrix(nrow=0,ncol=3)
geom1<-matrix(nrow=0,ncol=4)
dst1<-matrix(nrow=0,ncol=4)
all_previous_months<-0
timers<-c()
month_names<-c('January','February','March','April','May','June','July','August','September','October','November','December')
total_old<-0
total_new<-0
total_diff<-0
tec_dst_rejected<-matrix(nrow=0,ncol=2)
tec_dst_accepted<-matrix(nrow=0,ncol=2)
tec_dst_all<-matrix(nrow=0,ncol=2)
for(j in 1:12){ 
  month_file_name<-paste('darwin_data/data_all_',j,'.csv',sep='')
  data_month<-as.data.frame(read.csv(month_file_name,header=TRUE,sep=',',skip=0))
  old_len<-length(data_month$TEC)
  total_old<-total_old+old_len
  tec_dst_all_month<-cbind(data_month['TEC'],data_month['Dst'])
  colnames(tec_dst_all_month)<-c('TEC','Dst')
  tec_dst_all<-rbind(tec_dst_all,tec_dst_all_month)
  print(summary(tec_dst_all_month))
  data_rejected<-data_month[data_month$TEC >= limited,]
  data_month<-data_month[data_month$TEC < limited,]
  new_len<-length(data_month$TEC)
  total_new<-total_new+new_len
  diff_len<-old_len-new_len
  total_diff<-total_diff+diff_len
  tec_dst_rejected_month<-cbind(data_month['TEC'],data_month['Dst'])
  colnames(tec_dst_rejected_month)<-c('TEC','Dst')
  tec_dst_rejected<-rbind(tec_dst_rejected,tec_dst_rejected_month)
  print(summary(tec_dst_rejected_month))
  tec_dst_accepted_month<-cbind(data_rejected['TEC'],data_rejected['Dst'])
  colnames(tec_dst_accepted_month)<-c('TEC','Dst')
  tec_dst_accepted<-rbind(tec_dst_accepted,tec_dst_accepted_month)
  print(summary(tec_dst_accepted_month))
  print(paste(month_names[j],old_len,new_len,diff_len,round(diff_len/old_len*100,2),round(new_len/old_len*100,2)))
  day_value<-all_previous_months+data_month['d']
  timer<-day_value+data_month['total.seconds']/86400
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
  all_previous_months<-all_previous_months+length(unique(data_month[['d']]))
  print(max(diono$TEC))
}
print(paste(2014,total_old,total_new,total_diff,round(total_diff/total_old*100,2),round(total_new/total_old*100,2)))
print(max(iono1$TEC))
print(summary(tec_dst_rejected))
print(summary(tec_dst_accepted))
print(summary(tec_dst_all))

# Data input: Geomagnetic field density components,as taken at the Kakadu,NT
# INTERMAGNET reference station-original data reformatted

print(summary(geom1))
print(summary(iono1))
print(quantile(iono1$TEC,probs=seq(0,1,0.01)))
print(summary(dst1))

# TEC and B data aggregation into single data frame

envi<-merge(iono1,geom1,by='time')

### Model development

envi<-cbind(envi['TEC'],envi$Bx,envi$By,envi$Bz)
colnames(envi)<-c('TEC','Bx','By','Bz')
print(summary(envi))

# Define a 80%/20% train/test split of the data

trainIndex<-createDataPartition(envi$TEC,p=0.80,list=FALSE)
dataTrain<-envi[trainIndex,]
dataTest<-envi[-trainIndex,]
print(paste('envi', length(envi$TEC)))
print(paste('dataTrain', length(dataTrain$TEC)))
print(paste('dataTest', length(dataTest$TEC)))