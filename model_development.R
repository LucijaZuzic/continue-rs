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

# Stop writing to the file

sink()

# Start writing to an output file

limited<-150

dir.create(file.path(paste('TEC_',limited,sep='')),showWarnings=FALSE)
dir.create(file.path(paste('TEC_',limited,'/Klobuchar',sep='')),showWarnings=FALSE)
dir.create(file.path(paste('TEC_',limited,'/GAMB',sep='')),showWarnings=FALSE)
dir.create(file.path(paste('TEC_',limited,'/BCART',sep='')),showWarnings=FALSE)
dir.create(file.path(paste('TEC_',limited,'/CART1',sep='')),showWarnings=FALSE)
dir.create(file.path(paste('TEC_',limited,'/LM',sep='')),showWarnings=FALSE)
dir.create(file.path(paste('TEC_',limited,'/SGB',sep='')),showWarnings=FALSE)

sink(paste('TEC_',limited,'/analysis_output_',limited,'.txt',sep=''))

oldw<-getOption('warn')
options(warn=-1)

# Data input: TEC and dTEC data derived from RINEX GPS observations taken at Darwin,NT,Australia
# using Gopi Seemala's GPS TEC software - days in 2014

iono1<-matrix(nrow=0,ncol=3)
geom1<-matrix(nrow=0,ncol=4)
dst1<-matrix(nrow=0,ncol=4)
all_previous_months<-0
timers<-c()

for(j in 1:12){ 
  month_file_name<-paste('darwin_data/data_all_',j,'.csv',sep='')
  data_month<-as.data.frame(read.csv(month_file_name,header=TRUE,sep=',',skip=0))
  data_month<-data_month[data_month$TEC < limited,]
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

klobuchar_coeff<-as.data.frame(read.csv('Klobuchar.csv',header=TRUE,sep=',',skip=0))
iono1<-merge(iono1, klobuchar_coeff, by = "DOY", all.x = TRUE)

# Klobuchar model predictions development

ktec<-c()
month_ktec<-c()

c<-2.99792458E+08
RE<-6.378E+06
h<-350E+03
fip<-78.3 * pi / 180
lop<-291 * pi / 180
E<-90 * pi / 180 # zenit angle=0,to obtain VTEC
A<-0 * pi / 180 # heading north,irrelevant,since E=90
fiu<--12.8437111 * pi / 180 # latitude of Darwin,NT,Australia,-12.8437111 S
lou<-131.1327361 * pi / 180 # longitude of Darwin,NT,Australia,131.1327361 E

fi<-pi/2 - E - asin(RE*cos(E)/(RE + h))
fiI<-asin(sin(fiu)*cos(fi)+cos(fiu)*sin(fi)*cos(A))
loI<-lou + (fi*sin(A)/cos(fiI))
fim<-asin(sin(fiI)*sin(fip)+cos(fiI)*cos(fip)*cos(loI-lop))
use_plot<-FALSE

df_time_proper<-list()
df_month<-list()
for(j in 1:365){
  df_time_proper<-append(df_time_proper,list(iono1[iono1$DOY == j, ]$time_proper))
  df_month<-append(df_month,list(iono1[iono1$DOY == j, ]$month))
}

ktec_list<-vector("list", 365)
if (use_plot) {
  month_ktec_list<-vector("list", 365)
}

ptm<-Sys.time()
for(j in 1:365){
  alpha<-c(iono1$A1[j],iono1$A2[j],iono1$A3[j],iono1$A4[j])
  beta<-c(iono1$B1[j],iono1$B2[j],iono1$B3[j],iono1$B4[j])
  df_new<-df_time_proper[j][[1]]
  
  AI<-0
  for(k in 1:4){
    AI<-AI + alpha[k]*((fim/pi)^(k-1))
  }
  if(AI<0){AI<-0}
  
  PI<-0
  for(k in 1:4){
    PI<-PI + beta[k]*((fim/pi)^(k-1))
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
    plot(tr_vec,tec_vec,xlab='Time [days in 2014]',ylab='Klobuchar model equivalent ionospheric delay [m]',xlim=c(j,j + 1),main=paste('Klobuchar model equivalent ionospheric delay [m] for DOY',j,'in 2014'))
    dev.off()
  }
}
klobuchar.proc<-Sys.time() - ptm
print(paste('Processing time = ',as.numeric(klobuchar.proc),' s.',sep=''))

ktec<-do.call(rbind, ktec_list)
ktec<-as.data.frame(ktec)
ktec<-ktec$tec_vec
rtec<-iono1$TEC

if (use_plot) {
  month_ktec<-do.call(rbind, month_ktec_list)
  pdf(file=paste('TEC_',limited,'/TEC_time_series_plot_',limited,'.pdf',sep=''))
  it<-iono1$time
  plot(it,rtec,type='l',col='red',ylim=c(0,limited),ylab='TEC [TECU]',main='Klobuchar model and experimental (true) TEC in 2014',xlab='Days in 2014')
  lines(it,ktec,col='blue')
  legend('topright',inset=.005,c('experimental (true) TEC','Klobuchar model TEC'),fill=c('red','blue'),horiz=FALSE)
  dev.off()
  
  month_ktec<-as.data.frame(month_ktec)
  colnames(month_ktec)<-c('month','TEC')
  pdf(width=16,height=16,file=paste('TEC_',limited,'/TEC_month_time_series_plot_',limited,'.pdf',sep=''))
  par(mfrow=c(4,3),mar=c(5,5.5,4,1.5)+0.1)
  month_names<-c('January','February','March','April','May','June','July','August','September','October','November','December')
  for(j in 1:12){
    ktec2<-month_ktec[month_ktec$month == j,]$TEC
    rtec2<-iono1[iono1$month == j,]$TEC
    it2<-iono1[iono1$month == j,]$time
    plot(it2,rtec2,type='l',col='red',cex.main=2,cex.lab=2,cex.axis=2,ylim=c(0,limited),ylab='TEC [TECU]',main=paste('Klobuchar model and experimental (true) TEC\n',month_names[j],'2014'),xlab=paste('Days in',month_names[j],'2014'))
    lines(it2,ktec2,col='blue')
    if (j==12) {
      legend('topright',cex=2,inset=.005,c('experimental (true) TEC','Klobuchar model TEC'),fill=c('red','blue'),horiz=FALSE)
    }
  }
  par(mfrow=c(1,1))
  dev.off()
}

## Determine P-O, adjR2, RMSE, Q-Q for Klobuchar model

pdf(file=paste('TEC_',limited,'/Klobuchar/Klobuchar_plot_',limited,'.pdf',sep=''))
plot(rtec,ktec,main='Observed and Klobuchar model TEC [TECU]',xlab='Observed TEC [TECU]',ylab='Klobuchar model TEC [TECU]')
dev.off()
kl_testing<-as.numeric(postResample(pred=ktec,obs=rtec))
kl_MAE_testing<-kl_testing[3]
kl_adjR2_testing<-kl_testing[2]
kl_RMSE_testing<-kl_testing[1]
print(paste('RMSE', kl_RMSE_testing))
print(paste('adjR2', kl_adjR2_testing))
print(paste('MAE', kl_MAE_testing))
df_res<-data.frame(Residuals=ktec-rtec)
pdf(file = tempfile())
my_qq_k<-plot_qq(df_res)[[1]]+labs(title='Q-Q plot for Klobuchar model residuals',x='Theoretical normal quantiles',y='Sample quantiles (residuals)')
dev.off()
pdf(file=paste('TEC_',limited,'/Klobuchar/Klobuchar_QQ_plot_',limited,'.pdf',sep=''))
print(my_qq_k)
dev.off()
foo<-cbind(ktec,rtec)
colnames(foo)<-c('pred','real')
write.csv(foo,paste('TEC_',limited,'/Klobuchar/Klobuchar_',limited,'.csv',sep=''))

# Data input: Geomagnetic field density components,as taken at the Kakadu,NT
# INTERMAGNET reference station - original data reformatted

print(summary(geom1))

par(mfrow=c(1,1),mar=c(5,5.5,4,1.5)+0.1)
pdf(width=8,height=4.5,file=paste('TEC_',limited,'/Dst_plot_',limited,'.pdf',sep=''))
plot(dst1$time,dst1$Dst,type='l',col='red',main='Dst [nt] for days in 2014',xlab='Time [days in 2014]',ylab='Dst [nT]')
dev.off()

# TEC & B data aggregation into single data frame per geomagnetic event

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

# Prepare training scheme

trainControl<-trainControl(method='repeatedcv',number=10,repeats=5)

get_best_result=function(caret_fit) {
  best=which(rownames(caret_fit$results)==rownames(caret_fit$bestTune))
  best_result=caret_fit$results[best,]
  rownames(best_result)=NULL
  best_result
}

# Linear regression model

set.seed(13)
ptm<-Sys.time()
fit.lm<-lm(TEC~.,data=dataTrain)
lm.proc<-Sys.time() - ptm
print(paste('Processing time = ',as.numeric(lm.proc),' s.',sep=''))
print('Linear regression model')
print(summary(fit.lm))
lm_pred<-predict(fit.lm,dataTest)
lm.pr<-postResample(pred=lm_pred,obs=dataTest$TEC)
pdf(file=paste('TEC_',limited,'/LM/LM_plot_',limited,'.pdf',sep=''))
plot(dataTest$TEC,lm_pred,main='Observed and linear regression-predicted TEC [TECU]',xlab='Observed TEC [TECU]',ylab='Linear regression-predicted TEC [TECU]')
dev.off()
adjr2<-summary(fit.lm)$adj.r.squared
print(lm.pr)
df_res<-data.frame(Residuals=lm_pred-dataTest$TEC)
pdf(file = tempfile())
my_qq_lm<-plot_qq(df_res)[[1]]+labs(title='Q-Q plot for linear regression residuals',x='Theoretical normal quantiles',y='Sample quantiles (residuals)')
dev.off()
pdf(file=paste('TEC_',limited,'/LM/LM_QQ_plot_',limited,'.pdf',sep=''))
print(my_qq_lm)
dev.off()
foo<-cbind(lm_pred,dataTest$TEC)
colnames(foo)<-c('pred','real')
write.csv(foo,paste('TEC_',limited,'/LM/LM_',limited,'.csv',sep=''))

# CART1 

set.seed(13)
ptm<-Sys.time()
fit.cart1<-train(TEC~.,data=dataTrain,method='rpart1SE', trControl=trainControl)
cart1.proc<-Sys.time() - ptm
print(paste('Processing time = ',as.numeric(cart1.proc),' min.',sep=''))
cart1_pred<-predict(fit.cart1,dataTest)
cart1.pr<-postResample(pred=cart1_pred,obs=dataTest$TEC)
cart1<-get_best_result(fit.cart1)
pdf(file=paste('TEC_',limited,'/CART1/CART1_plot_',limited,'.pdf',sep=''))
plot(dataTest$TEC,cart1_pred,main='Observed and CART1-predicted TEC [TECU]',xlab='Observed TEC [TECU]',ylab='CART1-predicted TEC [TECU]')
dev.off()
print('CART1')
print(cart1.pr)
df_res<-data.frame(Residuals=cart1_pred-dataTest$TEC)
pdf(file = tempfile())
my_qq_cart1<-plot_qq(df_res)[[1]]+labs(title='Q-Q plot for CART1 residuals',x='Theoretical normal quantiles',y='Sample quantiles (residuals)')
dev.off()
pdf(file=paste('TEC_',limited,'/CART1/CART1_QQ_plot_',limited,'.pdf',sep=''))
print(my_qq_cart1)
dev.off()
foo<-cbind(cart1_pred,dataTest$TEC)
colnames(foo)<-c('pred','real')
write.csv(foo,paste('TEC_',limited,'/CART1/CART1_',limited,'.csv',sep=''))

# Boosted Generalized Additive Model

set.seed(13)
ptm<-Sys.time()
fit.gamb<-train(TEC~.,data=dataTrain,method='gamboost',trControl=trainControl)
gamb.proc<-Sys.time() - ptm
print(paste('Processing time = ',as.numeric(gamb.proc),' min.',sep=''))
gamb_pred<-predict(fit.gamb,dataTest)
gamb.pr<-postResample(pred=gamb_pred,obs=dataTest$TEC)
gamb<-get_best_result(fit.gamb)
pdf(file=paste('TEC_',limited,'/GAMB/GAMB_plot_',limited,'.pdf',sep=''))
plot(dataTest$TEC,gamb_pred,main='Observed and GAMB-predicted TEC [TECU]',xlab='Observed TEC [TECU]',ylab='GAMB-predicted TEC [TECU]')
dev.off()
print('Boosted Generalized Additive Model (GAMB)')
print(gamb.pr)
df_res<-data.frame(Residuals=gamb_pred-dataTest$TEC)
pdf(file = tempfile())
my_qq_gamb<-plot_qq(df_res)[[1]]+labs(title='Q-Q plot for GAMB residuals',x='Theoretical normal quantiles',y='Sample quantiles (residuals)')
dev.off()
pdf(file=paste('TEC_',limited,'/GAMB/GAMB_QQ_plot_',limited,'.pdf',sep=''))
print(my_qq_gamb)
dev.off()
foo<-cbind(gamb_pred,dataTest$TEC)
colnames(foo)<-c('pred','real')
write.csv(foo,paste('TEC_',limited,'/GAMB/GAMB_',limited,'.csv',sep=''))

## BAGGING

# Bagged CART

set.seed(13)
ptm<-Sys.time()
fit.bcart<-train(TEC~.,data=dataTrain,method='treebag',trControl=trainControl)
bcart.proc<-Sys.time() - ptm
print(paste('Processing time = ',as.numeric(bcart.proc),' min.',sep=''))
bcart_pred<-predict(fit.bcart,dataTest)
bcart.pr<-postResample(pred=bcart_pred,obs=dataTest$TEC)
bcart<-get_best_result(fit.bcart)
pdf(file=paste('TEC_',limited,'/BCART/BCART_plot_',limited,'.pdf',sep=''))
plot(dataTest$TEC,bcart_pred,main='Observed and BCART-predicted TEC [TECU]',xlab='Observed TEC [TECU]',ylab='BCART-predicted TEC [TECU]')
dev.off()
print('Bagged CART (BCART)')
print(bcart.pr)
df_res<-data.frame(Residuals=bcart_pred-dataTest$TEC)
pdf(file = tempfile())
my_qq_bcart<-plot_qq(df_res)[[1]]+labs(title='Q-Q plot for BCART residuals',x='Theoretical normal quantiles',y='Sample quantiles (residuals)')
dev.off()
pdf(file=paste('TEC_',limited,'/BCART/BCART_QQ_plot_',limited,'.pdf',sep=''))
print(my_qq_bcart)
dev.off()
foo<-cbind(bcart_pred,dataTest$TEC)
colnames(foo)<-c('pred','real')
write.csv(foo,paste('TEC_',limited,'/BCART/BCART_',limited,'.csv',sep=''))

# Stochastic Gradient Boosting

set.seed(13)
ptm<-Sys.time()
fit.sgb<-train(TEC~.,data=dataTrain,method='gbm',trControl=trainControl)
sgb.proc<-Sys.time() - ptm
print(paste('Processing time = ',as.numeric(sgb.proc),' min.',sep=''))
sgb_pred<-predict(fit.sgb,dataTest)
sgb.pr<-postResample(pred=sgb_pred,obs=dataTest$TEC)

sgb<-get_best_result(fit.sgb)
pdf(file=paste('TEC_',limited,'/SGB/SGB_plot_',limited,'.pdf',sep=''))
plot(dataTest$TEC,sgb_pred,main='Observed and SGB-predicted TEC [TECU]',xlab='Observed TEC [TECU]',ylab='SGB-predicted TEC [TECU]')
dev.off()
print('Stochastic Gradient Boosting (SGB)')
print(sgb.pr)
df_res<-data.frame(Residuals=sgb_pred-dataTest$TEC)
pdf(file = tempfile())
my_qq_sgb<-plot_qq(df_res)[[1]]+labs(title='Q-Q plot for SGB residuals',x='Theoretical normal quantiles',y='Sample quantiles (residuals)')
dev.off()
pdf(file=paste('TEC_',limited,'/SGB/SGB_QQ_plot_',limited,'.pdf',sep=''))
print(my_qq_sgb)
dev.off()
foo<-cbind(sgb_pred,dataTest$TEC)
colnames(foo)<-c('pred','real')
write.csv(foo,paste('TEC_',limited,'/SGB/SGB_',limited,'.csv',sep=''))

### Model performance assessment

par(mfrow=c(1,1))

# Create the data for the MAE chart

H<-c(lm.pr[3], cart1.pr[3],sgb.pr[3],bcart.pr[3],gamb.pr[3],kl_MAE_testing)
M<-c('LM','CART1','SGB','BCART','GAMB','K')

# Plot the bar chart

pdf(file=paste('TEC_',limited,'/MAE_',limited,'.pdf',sep=''))
bar_centers<-barplot(H,names.arg=M,xlab='Method',ylab='MAE [TECU]',ylim=c(0,max(H)*1.15),col='gray',main='MAE [TECU] for TEC predicted by each method',border='red')
text(x=bar_centers,y=H,labels=round(H,2),pos=3,cex=0.9,col="black")
dev.off()

par(mfrow=c(1,1))

# Create the data for the RMSE chart

H<-c(lm.pr[1],cart1.pr[1],sgb.pr[1],bcart.pr[1],gamb.pr[1],kl_RMSE_testing)
M<-c('LM','CART1','SGB','BCART','GAMB','K')

# Plot the bar chart

pdf(file=paste('TEC_',limited,'/RMSE_',limited,'.pdf',sep=''))
bar_centers<-barplot(H,names.arg=M,xlab='Method',ylab='RMSE [TECU]',ylim=c(0,max(H)*1.15),col='gray',main='RMSE [TECU] for TEC predicted by each method',border='red')
text(x=bar_centers,y=H,labels=round(H,2),pos=3,cex=0.9,col="black")
dev.off()

par(mfrow=c(1,1))

# Create the data for the adjR2 chart

H<-c(lm.pr[2]*100,cart1.pr[2]*100,sgb.pr[2]*100,bcart.pr[2]*100,gamb.pr[2]*100,kl_adjR2_testing*100)
M<-c('LM','CART1','SGB','BCART','GAMB','K')

# Plot the bar chart

pdf(file=paste('TEC_',limited,'/AdjR2_',limited,'.pdf',sep=''))
bar_centers<-barplot(H,names.arg=M,xlab='Method',ylab='Adjusted R2 [%]',ylim=c(0,max(H)*1.15),col='grey',main='Adjusted R2 [%] for TEC predicted by each method',border='red')
text(x=bar_centers,y=H,labels=round(H,2),pos=3,cex=0.9,col="black")
dev.off()

par(mfrow=c(1,1))

# Create the data for the model development time chart

H<-as.vector(as.numeric(c(lm.proc,cart1.proc,sgb.proc,bcart.proc,gamb.proc,klobuchar.proc)))
M<-c('LM','CART1','SGB','BCART','GAMB','K')

# Plot the bar chart

pdf(file=paste('TEC_',limited,'/Model_development_time_',limited,'.pdf',sep=''))
bar_centers<-barplot(H,names.arg=M,xlab='Method',ylab='Time elapsed [s]',ylim=c(0,max(H)*1.15),col='gray',main='Model development time [s]',border='red')
text(x=bar_centers,y=H,labels=round(H,2),pos=3,cex=0.9,col="black")
dev.off()

print(paste('MAE: LM = ',round(lm.pr[3],2),
            ', CART1 = ',round(cart1.pr[3],2),
            ', SGB = ',round(sgb.pr[3],2),
            ', BCART = ',round(bcart.pr[3],2),
            ', GAMB = ',round(gamb.pr[3],2),
            ', Klobuchar= ',round(kl_MAE_testing,2),sep=''))
print(paste('adj R2: LM = ',round(lm.pr[2]*100,2),
            '%, CART1 = ',round(cart1.pr[2]*100,2),
            '%, SGB = ',round(sgb.pr[2]*100,2),
            '%, BCART = ',round(bcart.pr[2]*100,2),
            '%, GAMB = ',round(gamb.pr[2]*100,2),
            '%, Klobuchar = ',round(kl_adjR2_testing*100,2),'%',sep=''))
print(paste('RMSE: LM = ',round(lm.pr[1],2),
            ', CART1 = ',round(cart1.pr[1],2),
            ', SGB = ',round(sgb.pr[1],2),
            ', BCART = ',round(bcart.pr[1],2),
            ', GAMB = ',round(gamb.pr[1],2),
            ', Klobuchar = ',round(as.numeric(kl_RMSE_testing),2),sep=''))
print(paste('development time: LM = ',round(as.numeric(lm.proc),2),' s',
            ', CART1 = ',round(as.numeric(cart1.proc),2),' min',
            ', SGB = ',round(as.numeric(sgb.proc),2),' min',
            ', BCART = ',round(as.numeric(bcart.proc),2),' min',
            ', GAMB = ',round(as.numeric(gamb.proc),2),' min',
            ', Klobuchar = ',round(as.numeric(klobuchar.proc),2),' s',sep=''))

pdf(file=paste('TEC_',limited,'/LM/LM_observed_predicted_',limited,'.pdf',sep=''))
par(mfrow=c(1,1))
plot(dataTest$TEC,lm_pred,main='Observed and linear regression-predicted TEC [TECU]',xlab='Observed TEC [TECU]',ylab='Linear regression-predicted TEC [TECU]')
a<-c(0,50)
b<-c(0,50)
dabl<-as.data.frame(cbind(a,b))
fit<-lm(b~a,data=dabl)
abline(fit,col='red')
dev.off()

pdf(file=paste('TEC_',limited,'/CART1/CART1_observed_predicted_',limited,'.pdf',sep=''))
par(mfrow=c(1,1))
plot(dataTest$TEC,cart1_pred,main='Observed and CART1-predicted TEC [TECU]',xlab='Observed TEC [TECU]',ylab='CART1-predicted TEC [TECU]')
a<-c(0,50)
b<-c(0,50)
dabl<-as.data.frame(cbind(a,b))
fit<-lm(b~a,data=dabl)
abline(fit,col='red')
dev.off()

pdf(file=paste('TEC_',limited,'/BCART/BCART_observed_predicted_',limited,'.pdf',sep=''))
par(mfrow=c(1,1))
plot(dataTest$TEC,bcart_pred,main='Observed and BCART-predicted TEC [TECU]',xlab='Observed TEC [TECU]',ylab='BCART-predicted TEC [TECU]')
a<-c(0,50)
b<-c(0,50)
dabl<-as.data.frame(cbind(a,b))
fit<-lm(b~a,data=dabl)
abline(fit,col='red')
dev.off()

pdf(file=paste('TEC_',limited,'/GAMB/GAMB_observed_predicted_',limited,'.pdf',sep=''))
par(mfrow=c(1,1))
plot(dataTest$TEC,gamb_pred,main='Observed and GAMB-predicted TEC [TECU]',xlab='Observed TEC [TECU]',ylab='GAMB-predicted TEC [TECU]')
a<-c(0,50)
b<-c(0,50)
dabl<-as.data.frame(cbind(a,b))
fit<-lm(b~a,data=dabl)
abline(fit,col='red')
dev.off()

pdf(file=paste('TEC_',limited,'/SGB/SGB_observed_predicted_',limited,'.pdf',sep=''))
par(mfrow=c(1,1))
plot(dataTest$TEC,sgb_pred,main='Observed and SGB-predicted TEC [TECU]',xlab='Observed TEC [TECU]',ylab='SGB-predicted TEC [TECU]')
a<-c(0,50)
b<-c(0,50)
dabl<-as.data.frame(cbind(a,b))
fit<-lm(b~a,data=dabl)
abline(fit,col='red')
dev.off()

pdf(file=paste('TEC_',limited,'/Klobuchar/Klobuchar_observed_predicted_',limited,'.pdf',sep=''))
par(mfrow=c(1,1))
plot(rtec,ktec,type='p',main='Observed and Klobuchar model TEC [TECU]',xlab='Observed TEC [TECU]',ylab='Klobuchar model TEC [TECU]')
a<-c(0,50)
b<-c(0,50)
dabl<-as.data.frame(cbind(a,b))
fit<-lm(b~a,data=dabl)
abline(fit,col='red')
dev.off()

pdf(file=paste('TEC_',limited,'/LM/LM_observed_predicted_',limited,'.pdf',sep=''))
par(mfrow=c(1,1))
plot(dataTest$TEC,lm_pred,main='Observed and linear regression-predicted TEC [TECU]',xlab='Observed TEC [TECU]',ylab='Linear regression-predicted TEC [TECU]')
a<-c(0,50)
b<-c(0,50)
dabl<-as.data.frame(cbind(a,b))
fit<-lm(b~a,data=dabl)
abline(fit,col='red')
dev.off()

pdf(width=16,height=9,file=paste('TEC_',limited,'/All_observed_predicted_',limited,'.pdf',sep=''))
par(mfrow=c(2,3))
plot(dataTest$TEC,lm_pred,main='Observed and linear regression-predicted TEC [TECU]',xlab='Observed TEC [TECU]',ylab='Linear regression-predicted TEC [TECU]')
a<-c(0,50)
b<-c(0,50)
dabl<-as.data.frame(cbind(a,b))
fit<-lm(b~a,data=dabl)
abline(fit,col='red')

plot(dataTest$TEC,cart1_pred,main='Observed and CART1-predicted TEC [TECU]',xlab='Observed TEC [TECU]',ylab='CART1-predicted TEC [TECU]')
a<-c(0,50)
b<-c(0,50)
dabl<-as.data.frame(cbind(a,b))
fit<-lm(b~a,data=dabl)
abline(fit,col='red')

plot(dataTest$TEC,sgb_pred,main='Observed and SGB-predicted TEC [TECU]',xlab='Observed TEC [TECU]',ylab='SGB-predicted TEC [TECU]')
a<-c(0,50)
b<-c(0,50)
dabl<-as.data.frame(cbind(a,b))
fit<-lm(b~a,data=dabl)
abline(fit,col='red')

plot(dataTest$TEC,bcart_pred,main='Observed and BCART-predicted TEC [TECU]',xlab='Observed TEC [TECU]',ylab='BCART-predicted TEC [TECU]')
a<-c(0,50)
b<-c(0,50)
dabl<-as.data.frame(cbind(a,b))
fit<-lm(b~a,data=dabl)
abline(fit,col='red')

plot(dataTest$TEC,gamb_pred,main='Observed and GAMB-predicted TEC [TECU]',xlab='Observed TEC [TECU]',ylab='GAMB-predicted TEC [TECU]')
a<-c(0,50)
b<-c(0,50)
dabl<-as.data.frame(cbind(a,b))
fit<-lm(b~a,data=dabl)
abline(fit,col='red')

plot(rtec,ktec,type='p',main='Observed and Klobuchar model TEC [TECU]',xlab='Observed TEC [TECU]',ylab='Klobuchar model TEC [TECU]')
a<-c(0,50)
b<-c(0,50)
dabl<-as.data.frame(cbind(a,b))
fit<-lm(b~a,data=dabl)
abline(fit,col='red')
par(mfrow=c(1,1))
dev.off()

pdf(width=16,height=9,file=paste('TEC_',limited,'/All_QQ_plot_',limited,'.pdf',sep=''))
par(mfrow=c(2,3))
grid.arrange(
  my_qq_lm, 
  my_qq_cart1, 
  my_qq_sgb, 
  my_qq_bcart, 
  my_qq_gamb, 
  my_qq_k, 
  nrow = 2, 
  ncol = 3
)
par(mfrow=c(1,1))
dev.off()

## fitdistrplus -> skewness,kurtosis diagram

pdf(file=paste('TEC_',limited,'/Klobuchar/Klobuchar_Cullen_Frey_plot_',limited,'.pdf',sep=''))
desc_k<-descdist(ktec,boot=1000)
mtext("Klobuchar model TEC [TECU]",side=3,line=0)
saved_plot<-recordPlot()
saveRDS(saved_plot,file=paste('TEC_',limited,'/Klobuchar/Klobuchar_Cullen_Frey_plot_',limited,'.rds',sep=''))
saveRDS(desc_k,file=paste('TEC_',limited,'/Klobuchar/Klobuchar_Cullen_Frey_plot_',limited,'_stats.rds',sep=''))
dev.off()

pdf(file=paste('TEC_',limited,'/Cullen_Frey_plot_',limited,'.pdf',sep=''))
desc_real<-descdist(iono1$TEC,boot=1000)
mtext("Observed TEC [TECU]",side=3,line=0)
saved_plot<-recordPlot()
saveRDS(saved_plot,file=paste('TEC_',limited,'/Cullen_Frey_plot_',limited,'.rds',sep=''))
saveRDS(desc_real,file=paste('TEC_',limited,'/Cullen_Frey_plot_',limited,'_stats.rds',sep=''))
dev.off()

pdf(file=paste('TEC_',limited,'/LM/LM_Cullen_Frey_plot_',limited,'.pdf',sep=''))
desc_lm<-descdist(lm_pred,boot=1000)
mtext("Linear-regression predicted TEC [TECU]",side=3,line=0)
saved_plot<-recordPlot()
saveRDS(saved_plot,file=paste('TEC_',limited,'/LM/LM_Cullen_Frey_plot_',limited,'.rds',sep=''))
saveRDS(desc_lm,file=paste('TEC_',limited,'/LM/LM_Cullen_Frey_plot_',limited,'_stats.rds',sep=''))
dev.off()

pdf(file=paste('TEC_',limited,'/CART1/CART1_Cullen_Frey_plot_',limited,'.pdf',sep=''))
desc_cart1<-descdist(cart1_pred,boot=1000)
mtext("CART1-predicted TEC [TECU]",side=3,line=0)
saved_plot<-recordPlot()
saveRDS(saved_plot,file=paste('TEC_',limited,'/CART1/CART1_Cullen_Frey_plot_',limited,'.rds',sep=''))
saveRDS(desc_cart1,file=paste('TEC_',limited,'/CART1/CART1_Cullen_Frey_plot_',limited,'_stats.rds',sep=''))
dev.off()

pdf(file=paste('TEC_',limited,'/GAMB/GAMB_Cullen_Frey_plot_',limited,'.pdf',sep=''))
desc_gamb<-descdist(as.numeric(gamb_pred),boot=1000)
mtext("GAMB-predicted TEC [TECU]",side=3,line=0)
saved_plot<-recordPlot()
saveRDS(saved_plot,file=paste('TEC_',limited,'/GAMB/GAMB_Cullen_Frey_plot_',limited,'.rds',sep=''))
saveRDS(desc_gamb,file=paste('TEC_',limited,'/GAMB/GAMB_Cullen_Frey_plot_',limited,'_stats.rds',sep=''))
dev.off()

pdf(file=paste('TEC_',limited,'/BCART/BCART_Cullen_Frey_plot_',limited,'.pdf',sep=''))
desc_bcart<-descdist(bcart_pred,boot=1000)
mtext("BCART-predicted TEC [TECU]",side=3,line=0)
saved_plot<-recordPlot()
saveRDS(saved_plot,file=paste('TEC_',limited,'/BCART/BCART_Cullen_Frey_plot_',limited,'.rds',sep=''))
saveRDS(desc_bcart,file=paste('TEC_',limited,'/BCART/BCART_Cullen_Frey_plot_',limited,'_stats.rds',sep=''))
dev.off()

pdf(file=paste('TEC_',limited,'/SGB/SGB_Cullen_Frey_plot_',limited,'.pdf',sep=''))
desc_sgb<-descdist(sgb_pred,boot=1000)
mtext("SGB-predicted TEC [TECU]",side=3,line=0)
saved_plot<-recordPlot()
saveRDS(saved_plot,file=paste('TEC_',limited,'/SGB/SGB_Cullen_Frey_plot_',limited,'.rds',sep=''))
saveRDS(desc_sgb,file=paste('TEC_',limited,'/SGB/SGB_Cullen_Frey_plot_',limited,'_stats.rds',sep=''))
dev.off()

pdf(width=16,height=9,file=paste('TEC_',limited,'/All_Cullen_Frey_plot_',limited,'.pdf',sep=''))
par(mfrow=c(2,3))
desc_lm<-descdist(lm_pred,boot=1000)
mtext("Linear-regression predicted TEC [TECU]",side=3,line=0)
desc_cart1<-descdist(cart1_pred,boot=1000)
mtext("CART1-predicted TEC [TECU]",side=3,line=0)
desc_sgb<-descdist(sgb_pred,boot=1000)
mtext("SGB-predicted TEC [TECU]",side=3,line=0)
desc_bcart<-descdist(bcart_pred,boot=1000)
mtext("BCART-predicted TEC [TECU]",side=3,line=0)
desc_gamb<-descdist(as.numeric(gamb_pred),boot=1000)
mtext("GAMB-predicted TEC [TECU]",side=3,line=0)
desc_k<-descdist(ktec,boot=1000)
mtext("Klobuchar model TEC [TECU]",side=3,line=0)
par(mfrow=c(1,1))
dev.off()

# Stop writing to the file
sink()
### END