rm(list=ls())

# Declaration of the R libraries utilised

# Core data frameworks
library(data.table)
library(tidyverse)

# Parallel processing
library(doParallel)

# Exploratory and diagnostic tools
library(DataExplorer)
library(fitdistrplus)

# Model wrapper
library(caret)

# Required by caret
library(glmnet)         # Generalized linear paths
library(relaimpo)       # Relative importance
library(monomvn)        # Parsimonious modeling
library(ranger)         # ranger for RF
library(gbm)            # gbm for SGB
library(nnet)           # nnet for NN

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

# Stop writing to the file

sink()
closeAllConnections()

# Start writing to an output file

limited<-150

# Synchronized arrays sorted logically from simplest baseline to heaviest ensemble

model_abbreviations<-c('LM','KNN','NN','RF','SGB')
model_names<-c('lm','knn','nnet','ranger','gbm')

dir.create(file.path(paste('TEC_',limited,sep='')),showWarnings=FALSE)
dir.create(file.path(paste('TEC_',limited,'/Klobuchar',sep='')),showWarnings=FALSE)

for (model_order in 1:length(model_abbreviations)) {
  dir.create(file.path(paste('TEC_',limited,'/',model_abbreviations[model_order],sep='')),showWarnings=FALSE)
}

sink(paste('TEC_',limited,'/analysis_output_',limited,'.txt',sep=''))

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

klobuchar_coeff<-as.data.frame(read.csv('Klobuchar.csv',header=TRUE,sep=',',skip=0))
iono1<-merge(iono1, klobuchar_coeff, by = "DOY", all.x = TRUE)

# Klobuchar model predictions development

ktec<-c()
month_ktec<-c()

c<-2.99792458E+08
RE<-6.378E+06
h<-350E+03
fip<-78.3*pi/180
lop<-291*pi/180
E<-90*pi/180 # zenit angle=0,to obtain VTEC
A<-0*pi/180 # heading north,irrelevant,since E=90
fiu<--12.8437111*pi/180 # latitude of Darwin,NT,Australia,-12.8437111 S
lou<-131.1327361*pi/180 # longitude of Darwin,NT,Australia,131.1327361 E

fi<-pi/2-E-asin(RE*cos(E)/(RE+h))
fiI<-asin(sin(fiu)*cos(fi)+cos(fiu)*sin(fi)*cos(A))
loI<-lou+(fi*sin(A)/cos(fiI))
fim<-asin(sin(fiI)*sin(fip)+cos(fiI)*cos(fip)*cos(loI-lop))
use_plot<-FALSE

df_time_proper<-list()
df_month<-list()
for(j in 1:365){
  df_time_proper<-append(df_time_proper,list(iono1[iono1$DOY==j, ]$time_proper))
  df_month<-append(df_month,list(iono1[iono1$DOY==j, ]$month))
}

ktec_list<-vector('list',365)
if (use_plot) {
  month_ktec_list<-vector('list',365)
}

ptm<-Sys.time()
for(j in 1:365){
  alpha<-c(iono1$A1[j],iono1$A2[j],iono1$A3[j],iono1$A4[j])
  beta<-c(iono1$B1[j],iono1$B2[j],iono1$B3[j],iono1$B4[j])
  df_new<-df_time_proper[j][[1]]
  
  AI<-0
  for(k in 1:4){
    AI<-AI+alpha[k]*((fim/pi)^(k-1))
  }
  if(AI<0){AI<-0}
  
  PI<-0
  for(k in 1:4){
    PI<-PI+beta[k]*((fim/pi)^(k-1))
  }
  if(PI<72000){PI<-72000}
  
  t_vec<-(43200*loI/pi+df_new) %% 86400
  
  XI_vec<-2*pi*(t_vec-50400)/PI
  
  I1_vec<-5E-09+(abs(XI_vec)<pi/2)*(AI*cos(XI_vec))
  
  tec_vec<-c*I1_vec*(1575.42e+06)^2/(40.3E+16)
  tr_vec<-j+df_new/86400
  
  ktec_list[[j]]<-cbind(tr_vec, tec_vec)
  if (use_plot) {
    df_month_new<-df_month[j][[1]]
    month_ktec_list[[j]]<-cbind(df_month_new,tec_vec)
    dir.create(file.path(paste('TEC_',limited,'/day_plots_',limited,sep='')),showWarnings=FALSE)
    pdf(file=paste('TEC_',limited,'/day_plots_',limited,'/Klobuchar_plot_',j,'_',limited,'.pdf',sep=''))
    plot(tr_vec,tec_vec,xlab='Time [days in 2014]',ylab='Klobuchar model equivalent ionospheric delay [m]',xlim=c(j,j+1),main=paste('Klobuchar model equivalent ionospheric delay [m] for DOY',j,'in 2014'))
    dev.off()
  }
}
klobuchar.proc<-as.numeric(difftime(Sys.time(),ptm,units='secs'))

ktec<-do.call(rbind, ktec_list)
ktec<-as.data.frame(ktec)
ktec<-ktec$tec_vec
rtec<-iono1$TEC

par(mfrow=c(1,1))
pdf(file=paste('TEC_',limited,'/TEC_Dst_plot_',limited,'.pdf',sep=''))
plot(iono1$TEC,dst1$Dst,main='TEC [TECU] and Dst [nT] for days in 2014',xlab='TEC [TECU]',ylab='Dst [nT]')
dev.off()

if (use_plot) {
  month_ktec<-do.call(rbind, month_ktec_list)
  pdf(file=paste('TEC_',limited,'/TEC_time_series_plot_',limited,'.pdf',sep=''))
  it<-iono1$time
  plot(it,rtec,type='l',col='red',ylim=c(0,max(rtec)),ylab='TEC [TECU]',main='Klobuchar model and experimental (true) TEC in 2014',xlab='Days in 2014')
  lines(it,ktec,col='blue')
  legend('topright',inset=.005,c('experimental (true) TEC','Klobuchar model TEC'),fill=c('red','blue'),horiz=FALSE)
  dev.off()
  
  month_ktec<-as.data.frame(month_ktec)
  colnames(month_ktec)<-c('month','TEC')
  pdf(width=16,height=16,file=paste('TEC_',limited,'/TEC_month_time_series_plot_',limited,'.pdf',sep=''))
  par(mfrow=c(4,3),mar=c(5,5.5,4,1.5)+0.1)
  month_names<-c('January','February','March','April','May','June','July','August','September','October','November','December')
  for(j in 1:12){
    ktec2<-month_ktec[month_ktec$month==j,]$TEC
    rtec2<-iono1[iono1$month==j,]$TEC
    it2<-iono1[iono1$month==j,]$time
    plot(it2,rtec2,type='l',col='red',cex.main=2,cex.lab=2,cex.axis=2,ylim=c(0,max(rtec2)),ylab='TEC [TECU]',main=paste('Klobuchar model and experimental (true) TEC\n',month_names[j],'2014'),xlab=paste('Days in',month_names[j],'2014'))
    lines(it2,ktec2,col='blue')
    if ((j==3&&limited==1290)||(j==7&&limited==150)) {
      legend('topright',cex=2,inset=.005,c('experimental (true) TEC','Klobuchar model TEC'),fill=c('red','blue'),horiz=FALSE)
    }
  }
  par(mfrow=c(1,1))
  dev.off()
}

## Determine P-O, adjR2, RMSE, Q-Q for Klobuchar model

pdf(file=paste('TEC_',limited,'/Klobuchar/Klobuchar_plot_',limited,'.pdf',sep=''))
plot(rtec,ktec,main='Observed and Klobuchar model TEC [TECU]',xlab='Observed TEC [TECU]',ylab='Klobuchar model TEC [TECU]')
a<-c(0,50)
b<-c(0,50)
dabl<-as.data.frame(cbind(a,b))
fit<-lm(b~a,data=dabl)
abline(fit,col='red')
dev.off()
kl_testing<-as.numeric(postResample(pred=ktec,obs=rtec))
kl_adjR2_testing<-1-(((1-kl_testing[2])*(length(ktec)-1))/(length(ktec)-1-1))
all_models_metrics<-list()
print('K')
print(paste('Testing time=',as.numeric(klobuchar.proc),' s.',sep=''))
print(kl_testing)
print(paste('Adjusted R-squared:',kl_adjR2_testing))

use_qq<-FALSE

if (use_qq) {
  df_res<-data.frame(Residuals=ktec-rtec)
  pdf(file=tempfile())
  my_qq_k<-plot_qq(df_res)[[1]]+labs(title='Q-Q plot for Klobuchar model residuals',x='Theoretical normal quantiles',y='Sample quantiles (residuals)')
  dev.off()
  pdf(file=paste('TEC_',limited,'/Klobuchar/Klobuchar_QQ_plot_',limited,'.pdf',sep=''))
  print(my_qq_k)
  dev.off()
}

foo<-cbind(ktec,rtec)
colnames(foo)<-c('pred','real')
write.csv(foo,paste('TEC_',limited,'/Klobuchar/Klobuchar_',limited,'.csv',sep=''))

## fitdistrplus -> skewness,kurtosis diagram

use_cullen<-FALSE

if (use_cullen) {
  pdf(file=paste('TEC_',limited,'/Cullen_Frey_plot_',limited,'.pdf',sep=''))
  desc_real<-descdist(iono1$TEC,boot=1000)
  mtext("Observed TEC [TECU]",side=3,line=0)
  saved_plot<-recordPlot()
  saveRDS(saved_plot,file=paste('TEC_',limited,'/Cullen_Frey_plot_',limited,'.rds',sep=''))
  saveRDS(desc_real,file=paste('TEC_',limited,'/Cullen_Frey_plot_',limited,'_stats.rds',sep=''))
  dev.off()
  
  pdf(file=paste('TEC_',limited,'/Klobuchar/Klobuchar_Cullen_Frey_plot_',limited,'.pdf',sep=''))
  desc_k<-descdist(ktec,boot=1000)
  mtext("Klobuchar model TEC [TECU]",side=3,line=0)
  saved_plot<-recordPlot()
  saveRDS(saved_plot,file=paste('TEC_',limited,'/Klobuchar/Klobuchar_Cullen_Frey_plot_',limited,'.rds',sep=''))
  saveRDS(desc_k,file=paste('TEC_',limited,'/Klobuchar/Klobuchar_Cullen_Frey_plot_',limited,'_stats.rds',sep=''))
  dev.off()
}

# Data input: Geomagnetic field density components,as taken at the Kakadu,NT
# INTERMAGNET reference station-original data reformatted

print(summary(geom1))
print(summary(iono1))
print(quantile(iono1$TEC,probs=seq(0,1,0.01)))
print(summary(dst1))

par(mfrow=c(1,1),mar=c(5,5.5,4,1.5)+0.1)
pdf(width=8,height=4.5,file=paste('TEC_',limited,'/Dst_plot_',limited,'.pdf',sep=''))
plot(dst1$time,dst1$Dst,type='l',col='red',main='Dst [nT] for days in 2014',xlab='Time [days in 2014]',ylab='Dst [nT]')
dev.off()

# TEC and B data aggregation into single data frame

envi<-merge(iono1,geom1,by='time')

pdf(file=paste('TEC_',limited,'/TEC_B_DOY_plot_',limited,'.pdf',sep=''))
par(mfrow=c(2,2))
plot(envi$time,envi$TEC,type='l',main='TEC [TECU] for days in 2014',col='red',xlab='Time [days in 2014]',ylab='TEC [TECU]')
plot(envi$time,envi$Bx,type='l',main='Bx [nT] for days in 2014',col='red',xlab='Time [days in 2014]',ylab='Bx [nT]')
plot(envi$time,envi$By,type='l',main='By [nT] for days in 2014',col='red',xlab='Time [days in 2014]',ylab='By [nT]')
plot(envi$time,envi$Bz,type='l',main='Bz [nT] for days in 2014',col='red',xlab='Time [days in 2014]',ylab='Bz [nT]')
par(mfrow=c(1,1))
dev.off()

pdf(file=paste('TEC_',limited,'/TEC_B_density_plot_',limited,'.pdf',sep=''))
par(mfrow=c(2,2))
plot(density(envi$TEC),main='Probability density TEC [TECU]',xlab='TEC [TECU]',ylab='Density')
plot(density(envi$Bx),main='Probability density Bx [nT]',xlab='Bx [nT]',ylab='Density')
plot(density(envi$By),main='Probability density By [nT]',xlab='By [nT]',ylab='Density')
plot(density(envi$Bz),main='Probability density Bz [nT]',xlab='Bz [nT]',ylab='Density')
par(mfrow=c(1,1))
dev.off()

pdf(file=paste('TEC_',limited,'/TEC_B_occurences_plot_',limited,'.pdf',sep=''))
par(mfrow=c(2,2))
boxplot(envi$TEC,main='Number of occurences TEC [TECU]',xlab='TEC [TECU]',ylab='Number of occurences')
boxplot(envi$Bx,main='Number of occurences Bx [nT]',xlab='Bx [nT]',ylab='Number of occurences')
boxplot(envi$By,main='Number of occurences By [nT]',xlab='By [nT]',ylab='Number of occurences')
boxplot(envi$Bz,main='Number of occurences Bz [nT]',xlab='Bz [nT]',ylab='Number of occurences')
par(mfrow=c(1,1))
dev.off()

pdf(file=paste('TEC_',limited,'/TEC_occurences_density_plot_',limited,'.pdf',sep=''))
par(mfrow=c(1,2))
boxplot(envi$TEC,main='Number of occurences TEC [TECU]',xlab='TEC [TECU]',ylab='Number of occurences')
plot(density(envi$TEC),main='Probability density TEC [TECU]',xlab='TEC [TECU]',ylab='Density')
par(mfrow=c(1,1))
dev.off()

pdf(file=paste('TEC_',limited,'/B_occurences_density_plot_',limited,'.pdf',sep=''))
par(mfrow=c(2,3))
boxplot(envi$Bx,main='Number of occurences Bx [nT]',xlab='Bx [nT]',ylab='Number of occurences')
boxplot(envi$By,main='Number of occurences By [nT]',xlab='By [nT]',ylab='Number of occurences')
boxplot(envi$Bz,main='Number of occurences Bz [nT]',xlab='Bz [nT]',ylab='Number of occurences')

plot(density(envi$Bx),main='Probability density Bx [nT]',xlab='Bx [nT]',ylab='Density')
plot(density(envi$By),main='Probability density By [nT]',xlab='By [nT]',ylab='Density')
plot(density(envi$Bz),main='Probability density Bz [nT]',xlab='Bz [nT]',ylab='Density')
dev.off()

### Model development

envi<-cbind(envi['TEC'],envi$Bx,envi$By,envi$Bz)
colnames(envi)<-c('TEC','Bx','By','Bz')
print(summary(envi))

# Define a 80%/20% train/test split of the data

set.seed(13)

trainIndex<-createDataPartition(envi$TEC,p=0.80,list=FALSE)
dataTrain<-envi[trainIndex,]
dataTest<-envi[-trainIndex,]

# Prepare training scheme

trainControl<-trainControl(method='repeatedcv',number=10,repeats=5)

cores<-detectCores()-1
cl<-makePSOCKcluster(cores)
registerDoParallel(cl)
for (model_order in 1:length(model_abbreviations)) {
  print(model_abbreviations[model_order])
  ptm<-Sys.time()
  if (model_names[model_order]=='ranger') {
    fit.caretmodel<-train(TEC~.,data=dataTrain,preProcess=c('scale','center'),method=model_names[model_order],trControl=trainControl(method='none'),tuneGrid=expand.grid(mtry=2,splitrule='variance',min.node.size=5),num.threads=1)
  } else {
    if (model_names[model_order]=='knn') {
      fit.caretmodel<-train(TEC~.,data=dataTrain,preProcess=c('scale','center'),method=model_names[model_order],trControl=trainControl(method='none'),tuneGrid=expand.grid(k=3))
    } else {
      if (model_names[model_order]=='nnet') {
        fit.caretmodel<-train(TEC~.,data=dataTrain,preProcess=c('scale','center'),method=model_names[model_order],trControl=trainControl(method='none'),tuneGrid=expand.grid(size=5,decay=0.01),linout=TRUE,maxit=300)
      } else {
        if (model_names[model_order]=='lm') {
          fit.caretmodel<-train(TEC~.,data=dataTrain,preProcess=c('scale','center'),method=model_names[model_order],trControl=trainControl(method='none'))
        } else {
          fit.caretmodel<-train(TEC~.,data=dataTrain,preProcess=c('scale','center'),method=model_names[model_order],trControl=trainControl)
        }
      }
    }
  }
  caretmodel.proc<-as.numeric(difftime(Sys.time(),ptm,units='secs'))
  print(paste('Training time=',as.numeric(caretmodel.proc),' s.',sep=''))
  ptm<-Sys.time()
  caretmodel_pred<-predict(fit.caretmodel,dataTest)
  caretmodel.pr<-postResample(pred=caretmodel_pred,obs=dataTest$TEC)
  caretmodel.test<-as.numeric(difftime(Sys.time(),ptm,units='secs'))
  print(paste('Testing time=',as.numeric(caretmodel.test),' s.',sep=''))
  print(caretmodel.pr)
  adj_r_squared<-1-(((1-caretmodel.pr[2])*(nrow(dataTest)-1))/(nrow(dataTest)-ncol(dataTest)-2))
  print(paste('Adjusted R-squared:',adj_r_squared))
  pdf(file=paste('TEC_',limited,'/',model_abbreviations[model_order],'/',model_abbreviations[model_order],'_plot_',limited,'.pdf',sep=''))
  plot(dataTest$TEC,caretmodel_pred,main=paste('Observed and ',model_abbreviations[model_order],'-predicted TEC [TECU]',sep=''),xlab='Observed TEC [TECU]',ylab=paste(model_abbreviations[model_order],'TEC [TECU]'))
  a<-c(0,50)
  b<-c(0,50)
  dabl<-as.data.frame(cbind(a,b))
  fit<-lm(b~a,data=dabl)
  abline(fit,col='red')
  dev.off()
  if (use_qq) {
    df_res<-data.frame(Residuals=caretmodel_pred-dataTest$TEC)
    pdf(file=tempfile())
    my_qq_caretmodel<-plot_qq(df_res)[[1]]+labs(title=paste('Q-Q plot for',model_abbreviations[model_order],'residuals'),x='Theoretical normal quantiles',y='Sample quantiles (residuals)')
    dev.off()
    pdf(file=paste('TEC_',limited,'/',model_abbreviations[model_order],'/',model_abbreviations[model_order],'_QQ_plot_',limited,'.pdf',sep=''))
    print(my_qq_caretmodel)
    dev.off()
  }
  foo<-cbind(caretmodel_pred,dataTest$TEC)
  colnames(foo)<-c('pred','real')
  write.csv(foo,paste('TEC_',limited,'/',model_abbreviations[model_order],'/',model_abbreviations[model_order],'_',limited,'.csv',sep=''))
  if (use_cullen) {
    pdf(file=paste('TEC_',limited,'/',model_abbreviations[model_order],'/',model_abbreviations[model_order],'_Cullen_Frey_plot_',limited,'.pdf',sep=''))
    desc_caretmodel<-descdist(as.numeric(caretmodel_pred),boot=1000)
    mtext(paste(model_abbreviations[model_order],'-predicted TEC [TECU]',sep=''),side=3,line=0)
    saved_plot<-recordPlot()
    saveRDS(saved_plot,file=paste('TEC_',limited,'/',model_abbreviations[model_order],'/',model_abbreviations[model_order],'_Cullen_Frey_plot_',limited,'.rds',sep=''))
    saveRDS(desc_caretmodel,file=paste('TEC_',limited,'/',model_abbreviations[model_order],'/',model_abbreviations[model_order],'_Cullen_Frey_plot_',limited,'_stats.rds',sep=''))
    dev.off()
  }
  current_model_stats<-as.data.frame(t(caretmodel.pr))
  colnames(current_model_stats)<-c('RMSE','Rsquared','MAE')
  current_model_stats$training<-caretmodel.proc
  current_model_stats$testing<-caretmodel.test
  current_model_stats$adjRsquared<-adj_r_squared
  current_model_stats$name<-model_abbreviations[model_order]
  all_models_metrics[[model_order]]<-current_model_stats
  final_metrics_df<-do.call(rbind,all_models_metrics)
  final_metrics_df<-final_metrics_df[,c('name','RMSE','Rsquared','adjRsquared','MAE','training','testing')]
  master_csv_path<-paste('TEC_',limited,'/All_Models_Metrics_',limited,'.csv',sep='')
  write.csv(final_metrics_df, master_csv_path,row.names=FALSE)
  print(paste('Saved master metrics to:', master_csv_path))
}
stopCluster(cl)
registerDoSEQ()

current_model_stats<-as.data.frame(t(kl_testing))
colnames(current_model_stats)<-c('RMSE','Rsquared','MAE')
current_model_stats$training<-NA
current_model_stats$testing<-klobuchar.proc
current_model_stats$adjRsquared<-kl_adjR2_testing
current_model_stats$name<-'K'
all_models_metrics[[length(model_abbreviations)+1]]<-current_model_stats
final_metrics_df<-do.call(rbind,all_models_metrics)
final_metrics_df<-final_metrics_df[,c('name','RMSE','Rsquared','adjRsquared','MAE','training','testing')]
master_csv_path<-paste('TEC_',limited,'/All_Models_Metrics_',limited,'.csv',sep='')
write.csv(final_metrics_df, master_csv_path,row.names=FALSE)
print(paste('Saved master metrics to:', master_csv_path))

# Stop writing to the file

sink()
closeAllConnections()