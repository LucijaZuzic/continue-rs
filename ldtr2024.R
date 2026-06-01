rm(list=ls())

#Set path to the working directory containing:
#1. INTERMAGNET data set - to be taken from the INTERMAGNET web-site, and reformatted 
#accordingly - #reformatted original data enclosed
#2. TEC data set, as derived from the RINEX GPS observations using 
# GPS TEC software - enclosed

setwd('E:/Academic/0 Research/000 In preparation/000 A Dst-based space weather conditions machine learning classification model for GNSS PNT performance analysis/')

#Declaration of the R libraries utilised

library(data.table)
library(DataExplorer)
library(corrplot)
library(lubridate)
library(fitdistrplus)
library(weird)

oldw <- getOption("warn")
options(warn = -1)

# Reading the annual observed TEC data at Darwin, NT in 2014 and
# data frame arrangement

iono <- matrix(nrow=0, ncol = 3)
data <- as.data.frame(read.csv('darwin2014.csv', header=TRUE, 
                               sep =',', skip=0))

data <- data[,-1]
print(summary(data))
pe <- data[data$min==0 & data$sec==0,]

# Reading the annual estimated Dst and ap data in 2014 and
# data frame arrangement

geo <- as.data.frame(read.csv('Dstap2014.dat',header=FALSE,sep=''))
colnames(geo) <- c('year','DOY','hr','Dst','ap')
summary(geo)

# Merging TEC and Dst&ap data frames
dataset <- merge(pe, geo, by = c('DOY', 'hr'))
dataset <- dataset[,-c(3,4,10)]
summary(dataset)

# Exloratory analysis of observations per variables
plot_histogram(dataset, title = 'over-all TEC')

# Data set reduction by expelling TEC>=300 outliers
iono2 <- dataset[dataset$TEC<300,]
iono2 <- dataset[c(1:8000),]
plot_histogram(iono2, title = 'TEC<300')
#boxplot(dataset, col = 'red', names = c('overall TEC','TEC<200'), 
#ylab='values [TECU]', cex.axis=1.5, cex.lab=1.5)
dataset2 <- dataset[,-c(1,2)]
plot_correlation(dataset2, maxcat = 5L)
#Box-plots per Dst values, source of the classification criteria
plot_boxplot(dataset2, by = "Dst")
iono3 <- iono2[,-c(1,2)]
plot_boxplot(iono3, by = "Dst")
plot_scatterplot(dataset2, by = 'TEC', sampled_rows = 1000L)
plot_scatterplot(iono3, by = 'Dst', sampled_rows = 1000L)

### Dst-based classification under the following rule:
# 15 ... 50 positive phase of the storm -> P
# -20 ... 15 normal -> N
# -55 ... -20 recovery phase (mostly) -> R
# -85 ... -55 through -> T
# -120 ... -85 extreme -> E


Dst_class <- character()
P <- iono3[iono3$Dst >= 15,]
N <- iono3[(iono3$Dst < 15 & iono3$Dst >= -20),]
R <- iono3[(iono3$Dst < (-20) & iono3$Dst >= (-55)),]
T <- iono3[(iono3$Dst < (-55) & iono3$Dst >= (-85)),]
E <- iono3[(iono3$Dst <= -85),]

# Adding estimate class (a categorical variable) to the data set
for(i in 1:length(iono3$TEC)){
  if(iono3$Dst[i]>=15){
    Dst_class[i] <- 'P'
  }else if(iono3$Dst[i]< 15 & iono3$Dst[i]>= -20){
    Dst_class[i] <- 'N'
  }else if(iono3$Dst[i]< -20 & iono3$Dst[i]>= -55){
    Dst_class[i] <- 'R'
  }else if(iono3$Dst[i]< -55 & iono3$Dst[i]>= -85){
    Dst_class[i] <- 'T'
  }else if(iono3$Dst[i]<= -85){
    Dst_class[i] <- 'E'
  }
}
classdata <- as.data.frame(cbind(iono3,Dst_class))


### Classification model development

library(caret)
library(mlbench)
library(MLmetrics)
#data(Sonar)

# Data set splitting into training and testing set
set.seed(13)
inTrain <- createDataPartition(
  y = classdata$Dst_class,
  ## the outcome data are needed
  p = .8, # partition value: .8 = 80%
  ## The percentage of data in the
  ## training set
  list = FALSE
)

str(inTrain)

training <- classdata[ inTrain,]
testing  <- classdata[-inTrain,]

nrow(training)
nrow(testing)

ctrl <- trainControl(
  method = "repeatedcv",
  number = 10,
  repeats = 10,
  classProbs = TRUE, 
  summaryFunction = multiClassSummary
)

## Development of candidate models 

# Candidate  model 1: svmPoly
# Least Squares Support Vector Machine with Polynomial Kernel
# Source: https://topepo.github.io/caret/train-models-by-tag.html

SVModel <- train(Dst_class ~ ., data = training,
                 method = "svmPoly",
                 trControl= ctrl,
                 tuneGrid = data.frame(degree = 1,
                                       scale = 1,
                                       C = 1),
                 preProcess = c("pca","scale","center"),
                 na.action = na.omit
)

SVMPredictions <-predict(SVModel, testing)
# Create confusion matrix
trueclasses <- as.factor(testing$Dst_class)
cmSVM <-confusionMatrix(SVMPredictions, trueclasses)
print(cmSVM)

importance <- varImp(SVModel, scale=FALSE)
plot(importance)


# Candidate  model 2: C5.0
# Decision Tree C5.0
# Source: https://topepo.github.io/caret/train-models-by-tag.html

DecTreeModel <- train(Dst_class ~ ., data = training, 
                      method = "C5.0",
                      preProcess=c("scale","center"),
                      trControl= ctrl,
                      na.action = na.omit
)

#Predictions
DTPredictions <-predict(DecTreeModel, testing, na.action = na.pass)
# Print confusion matrix and results
trueclasses <- as.factor(testing$Dst_class)
cmTree <-confusionMatrix(DTPredictions, trueclasses)
print(cmTree)


# Candidate  model 3: nb
# Naive Bayes
# Source: https://topepo.github.io/caret/train-models-by-tag.html
NaiveModel <- train(training, training$Dst_class, 
                    method = "nb",
                    preProcess=c("scale","center"),
                    trControl= ctrl,
                    na.action = na.omit
)

#Predictions
NaivePredictions <-predict(NaiveModel, testing, na.action = na.pass)
trueclasses <- as.factor(testing$Dst_class)
cmNaive <-confusionMatrix(NaivePredictions, trueclasses)
print(cmNaive)

# train model with neural networks
NNModel <- train(training, training$Dst_class,
                 method = "nnet",
                 trControl= ctrl,
                 preProcess=c("scale","center"),
                 na.action = na.omit
)

NNPredictions <-predict(NNModel, testing)
# Create confusion matrix
trueclasses <- as.factor(testing$Dst_class)
cmNN <-confusionMatrix(NNPredictions, trueclasses)
print(cmNN)



# Create four models as an ensemble
#library(caretEnsemble)
#library(caret)
#econtrol <- trainControl(method="cv", number=10, 
#                         savePredictions=TRUE, classProbs=TRUE)
#model_list <- caretList(Dst_class ~., data=training,
#                        methodList=c("svmPoly", "nnet", "C5.0", "nb"),
#                        preProcess=c("scale","center"),
#                        trControl = econtrol
#)

#results <- resamples(model_list)
# What is model correlation?
#mcr <-modelCor(results)
#print (mcr)

# Candidate  model 3: nb
# Naive Bayes
# Source: https://topepo.github.io/caret/train-models-by-tag.html

# Candidate  model 4: pls
# Partial Least Squares
# Source: https://topepo.github.io/caret/train-models-by-tag.html

plsFit <- train(
  Dst_class ~ .,
  data = training,
  method = "pls",
  preProc = c("center", "scale"),
  tuneLength = 15,
  trControl = ctrl,
  metric = "ROC"
)
plsFit

ggplot(plsFit)
plsClasses <- predict(plsFit, newdata = testing)
str(plsClasses)

plsProbs <- predict(plsFit, newdata = testing, type = "prob")
head(plsProbs)

trueclasses <- as.factor(testing$Dst_class)
confusionMatrix(data = plsClasses, trueclasses)

# Candidate  model 5: fda
# Flexible Discriminant Analysis
# Source: https://topepo.github.io/caret/train-models-by-tag.html
library(earth,mda)
fdaModel <- train(Dst_class ~ ., data = training,
                 method = "fda",
                 trControl= ctrl,
                 preProcess = c("pca","scale","center"),
                 na.action = na.omit
)

fdaPredictions <-predict(fdaModel, testing)
# Create confusion matrix
trueclasses <- as.factor(testing$Dst_class)
cmfda <-confusionMatrix(fdaPredictions, trueclasses)
print(cmfda)

importance <- varImp(SVModel, scale=FALSE)
plot(importance)


# Candidate  model 6: glmStepAIC
# Flexible Discriminant Analysis
# Source: https://topepo.github.io/caret/train-models-by-tag.html
library(nnet)
fdaModel <- train(Dst_class ~ ., data = training,
                  method = "pcaNNet",
                  trControl= ctrl,
                  preProcess = c("scale","center"),
                  na.action = na.omit
)

fdaPredictions <-predict(fdaModel, testing)
# Create confusion matrix
trueclasses <- as.factor(testing$Dst_class)
cmfda <-confusionMatrix(fdaPredictions, trueclasses)
print(cmfda)

importance <- varImp(fdaModel, scale=FALSE)
plot(importance)

# Candidate  model 6: glmStepAIC
# Flexible Discriminant Analysis
# Source: https://topepo.github.io/caret/train-models-by-tag.html

library(binda)

fdaModel <- train(Dst_class ~ ., data = training,
                  method = "binda",
                  trControl= ctrl,
                  preProcess = c("scale","center"),
                  na.action = na.omit
)

fdaPredictions <-predict(fdaModel, testing)
# Create confusion matrix
trueclasses <- as.factor(testing$Dst_class)
cmfda <-confusionMatrix(fdaPredictions, trueclasses)
print(cmfda)

importance <- varImp(fdaModel, scale=FALSE)
plot(importance)

#####################################################