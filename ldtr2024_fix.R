rm(list=ls())

#Set path to the working directory containing:
#1. INTERMAGNET data set - to be taken from the INTERMAGNET web-site, and reformatted 
#accordingly - #reformatted original data enclosed
#2. TEC data set, as derived from the RINEX GPS observations using 
# GPS TEC software - enclosed

library(tidyverse)
getCurrentFileLocation <-  function()
{
    this_file <- commandArgs() %>% 
    tibble::enframe(name = NULL) %>%
    tidyr::separate(col=value, into=c("key", "value"), sep="=", fill='right') %>%
    dplyr::filter(key == "--file") %>%
    dplyr::pull(value)
    if (length(this_file)==0)
    {
      this_file <- rstudioapi::getSourceEditorContext()$path
    }
    return(dirname(this_file))
}

setwd(getCurrentFileLocation())

#Declaration of the R libraries utilised

library(data.table)
library(DataExplorer)
library(corrplot)
library(lubridate)
library(fitdistrplus)
library(weird)

oldw <- getOption("warn")
options(warn = -1)

sink("data_vars.txt")

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
dev.copy(png, filename = "allTEC.png")
dev.off()
# Data set reduction by expelling TEC>=300 outliers
iono2 <- dataset[dataset$TEC<300,]
iono2 <- dataset[c(1:8000),]
plot_histogram(iono2, title = 'TEC<300')
dev.copy(png, filename = "300TEC.png")
dev.off()
#boxplot(dataset, col = 'red', names = c('overall TEC','TEC<200'), 
#ylab='values [TECU]', cex.axis=1.5, cex.lab=1.5)
dataset2 <- dataset[,-c(1,2)]
plot_correlation(dataset2, maxcat = 5L)
dev.copy(png, filename = "correlation.png")
dev.off()
#Box-plots per Dst values, source of the classification criteria
plot_boxplot(dataset2, by = "Dst")
dev.copy(png, filename = "dataset2boxplot.png")
dev.off()
iono3 <- iono2[,-c(1,2)]
plot_boxplot(iono3, by = "Dst")
dev.copy(png, filename = "iono3boxplot.png")
dev.off()
plot_scatterplot(dataset2, by = 'TEC', sampled_rows = 1000L)
dev.copy(png, filename = "dataset2scatterplot.png")
dev.off()
plot_scatterplot(iono3, by = 'Dst', sampled_rows = 1000L)
dev.copy(png, filename = "iono3scatterplot.png")
dev.off()

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
classdata <- as.data.frame(cbind(iono3, Dst_class))

classdata_no_Dst <- classdata[c("ap", "Bx", "By", "Bz", "TEC", "dTEC", "Dst_class")]
colnames(classdata_no_Dst) <- c("ap", "Bx", "By", "Bz", "TEC", "dTEC", "Dst_class")

classdata_no_TEC <- classdata[c("ap", "Bx", "By", "Bz", "Dst_class")]
colnames(classdata_no_TEC) <- c("ap", "Bx", "By", "Bz", "Dst_class")

classdata_coord <- classdata[c("Bx", "By", "Bz", "Dst_class")]
colnames(classdata_coord) <- c("Bx", "By", "Bz", "Dst_class")

classdata_xy_ap <- classdata[c("ap", "Bx", "By", "Dst_class")]
colnames(classdata_xy_ap) <- c("ap", "Bx", "By", "Dst_class")

classdata_xz_ap <- classdata[c("ap", "Bx", "Bz", "Dst_class")]
colnames(classdata_xz_ap) <- c("ap", "Bx", "Bz", "Dst_class")

classdata_yz_ap <- classdata[c("ap", "By", "Bz", "Dst_class")]
colnames(classdata_yz_ap) <- c("ap", "By", "Bz", "Dst_class")

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

training_all <- classdata[ inTrain,]
testing_all  <- classdata[-inTrain,]

training_no_Dst <- classdata_no_Dst[ inTrain,]
testing_no_Dst <- classdata_no_Dst[-inTrain,]
  
training_no_TEC <- classdata_no_TEC[ inTrain,]
testing_no_TEC <- classdata_no_TEC[-inTrain,]
  
training_coord <- classdata_coord[ inTrain,]
testing_coord <- classdata_coord[-inTrain,]

training_xy_ap <- classdata_xy_ap[ inTrain,]
testing_xy_ap <- classdata_xy_ap[-inTrain,]

training_xz_ap <- classdata_xz_ap[ inTrain,]
testing_xz_ap <- classdata_xz_ap[-inTrain,]

training_yz_ap <- classdata_yz_ap[ inTrain,]
testing_yz_ap <- classdata_yz_ap[-inTrain,]

nrow(training_all)
nrow(testing_all)

ctrl <- trainControl(
  method = "repeatedcv",
  number = 10,
  repeats = 10,
  classProbs = TRUE, 
  summaryFunction = multiClassSummary
)
sink()

## Development of candidate models 
# Source: https://topepo.github.io/caret/train-models-by-tag.html

library(earth,mda) # only for fda
library(nnet) # only for pcaNNet
library(binda) # only for binda

# Candidate model 1: svmPoly
# Least Squares Support Vector Machine with Polynomial Kernel
# Candidate model 2: C5.0
# Decision Tree C5.0
# Candidate model 3: nb
# Naive Bayes
# Candidate model 4: nnet
# Neural Networks
# Candidate  model 5: pls
# Partial Least Squares
# Candidate model 6: fda
# Flexible Discriminant Analysis
# Candidate model 7: pcaNNet
# Candidate model 8: binda
# Candidate model 9: glmStepAIC ?

#model_name_list <- list("svmPoly", "C5.0", "nb", "nnet", "pls", "fda", "pcaNNet", "binda", "glmStepAIC")
model_name_list <- list("svmPoly", "C5.0", "nb", "nnet", "pls", "fda", "pcaNNet")

model_trained_global <- 0
model_predictions_global <- 0
testing <- 0
training <- 0
trueclasses <- 0
  
my_train <- function(model_name) {
  if (model_name == "svmPoly") {
    model_trained <- train(Dst_class ~ ., data = training,
                           method = "svmPoly",
                           trControl= ctrl,
                           tuneGrid = data.frame(degree = 1,
                                                 scale = 1,
                                                 C = 1),
                           preProcess = c("pca","scale","center"),
                           na.action = na.omit
    )
  } else if (model_name == "C5.0") {
    model_trained <- train(Dst_class ~ ., data = training, 
                           method = "C5.0",
                           preProcess=c("scale","center"),
                           trControl= ctrl,
                           na.action = na.omit
    )
  } else if (model_name == "nb") {
    model_trained <- train(training, training$Dst_class, 
                           method = "nb",
                           preProcess=c("scale","center"),
                           trControl= ctrl,
                           na.action = na.omit
    )
  } else if (model_name == "nnet") {
    model_trained <- train(training, training$Dst_class,
                           method = "nnet",
                           trControl= ctrl,
                           preProcess=c("scale","center"),
                           na.action = na.omit
    )
  } else if (model_name == "pls") {
    model_trained <- train(
      Dst_class ~ .,
      data = training,
      method = "pls",
      preProc = c("center", "scale"),
      tuneLength = 15,
      trControl = ctrl,
      metric = "ROC"
    )
    print(model_trained)
  } else if (model_name == "fda") {
    model_trained <- train(Dst_class ~ ., data = training,
                           method = "fda",
                           trControl= ctrl,
                           preProcess = c("pca","scale","center"),
                           na.action = na.omit
    )
  } else if (model_name == "pcaNNet") {
    model_trained <- train(Dst_class ~ ., data = training,
                           method = "pcaNNet",
                           trControl= ctrl,
                           preProcess = c("scale","center"),
                           na.action = na.omit
    )
  } else if (model_name == "binda") {
    model_trained  <- train(Dst_class ~ ., data = training,
                            method = "binda",
                            trControl= ctrl,
                            preProcess = c("scale","center"),
                            na.action = na.omit
    )
  } else if (model_name == "glmStepAIC") {
    model_trained  <- train(Dst_class ~ ., data = training,
                            method = "glmStepAIC",
                            trControl= ctrl,
                            preProcess = c("scale","center"),
                            na.action = na.omit
    )
  }
  assign("model_trained_global", model_trained, envir = .GlobalEnv)
}

my_predict <- function(model_name_predict, model_type_predict) {
  # Predictions
  if (model_name_predict == "pls") {
    ggplot(model_trained_global)
    ggsave(paste(model_type_predict, "pls_plot.png", sep ="/"))
    model_predictions <- predict(model_trained_global, newdata = testing)
    print(str(model_predictions))
    model_probs <- predict(model_trained_global, newdata = testing, type = "prob")
    print(head(model_probs))
  } 
  model_predictions <- predict(model_trained_global, testing)
  # Create confusion matrix
  cm_model <- confusionMatrix(model_predictions, trueclasses)
  # Print confusion matrix and results
  print(cm_model)
  assign("model_predictions_global", model_predictions, envir = .GlobalEnv)
}

my_importance <- function(model_name_importance, model_type_importance) {
  # Variable importance
  importance <- varImp(model_trained_global, scale = FALSE)
  print(importance)
  # Spremanje dijagrama
  #plot(importance)
  #dev.copy(png, filename = paste(paste(paste(model_type_importance, "importance", sep = "/"), model_name_importance, sep = "_"), "png", sep = "."))
  #dev.off()
}
  
my_process <- function(model_name, model_type) {
  time_taken_train <- system.time(my_train(model_name))
  print("Train time")
  print(time_taken_train)
  time_taken_predict <- system.time(my_predict(model_name, model_type))
  print("Predict time")
  print(time_taken_predict)
  if (model_name == "pls" || model_name == "fda" || model_name == "C5.0") {
    time_taken_importance <- system.time(my_importance(model_name, model_type))
    print("Importance time")
    print(time_taken_importance)
  }
}

# Create four models as an ensemble
library(caretEnsemble)

my_ansamble <- function(list_all) {
  econtrol <- trainControl(method="cv", number=10, 
                           savePredictions=TRUE, classProbs=TRUE)
  model_list <- caretList(Dst_class ~., data=training,
                          methodList=list_all,
                          preProcess=c("scale","center"),
                          trControl = econtrol
  )
  results <- resamples(model_list)
  # What is model correlation?
  mcr <- modelCor(results)
  print(mcr)
}

model_types <- c("all", "no_Dst", "no_TEC", "coord", "xyap", "xzap", "yzap")

for (i in 1:length(model_types)) {
  model_type_use <- model_types[[i]]
  if (model_type_use == "all") {
    testing <- testing_all
    training <- training_all
  } else if (model_type_use == "no_Dst") {
    testing <- testing_no_Dst
    training <- training_no_Dst
  } else if (model_type_use == "no_TEC") {
    testing <- testing_no_TEC
    training <- training_no_TEC
  } else if (model_type_use == "coord") {
    testing <- testing_coord
    training <- training_coord
  } else if (model_type_use == "xyap") {
    testing <- testing_xy_ap
    training <- training_xy_ap
  } else if (model_type_use == "xzap") {
    testing <- testing_xz_ap
    training <- training_xz_ap
  } else if (model_type_use == "yzap") {
    testing <- testing_yz_ap
    training <- training_yz_ap
  }
  if (!file.exists(model_type_use)){
    dir.create(file.path(getCurrentFileLocation(), model_type_use))
  }
  trueclasses <- as.factor(testing$Dst_class)
  for (j in 1:length(model_name_list)) {
    model_name_use <- model_name_list[[j]]
    file_name <- paste(paste(model_type_use, model_name_use, sep = "/"), "txt", sep = ".")
    if (file.exists(file_name)) {
      print(model_type_use)
      print(model_name_use)
      sink(file_name)
      print(model_type_use)
      print(model_name_use)
      time_taken <- system.time(my_process(model_name_use, model_type_use))
      print("Total time")
      print(time_taken)
      sink()
      print(time_taken)
    }
  }
  
  list4 <- c("svmPoly", "nnet", "C5.0", "nb")
  list7 <- c("svmPoly", "C5.0", "nb", "nnet", "pls", "fda", "pcaNNet")
  
  if (file.exists(paste(model_type_use, "ansamble4.txt", sep = "/"))) {
      print(model_type_use)
      print("ansamble4")
      sink(paste(model_type_use, "ansamble4.txt", sep = "/"))
      print(model_type_use)
      print("ansamble4")
      time_taken <- system.time(my_ansamble(list4))
      print(time_taken)
      sink()
      print(time_taken)
  }
  
  if (file.exists(paste(model_type_use, "ansamble7.txt", sep = "/"))) {
      print(model_type_use)
      print("ansamble7")
      sink(paste(model_type_use, "ansamble7.txt", sep = "/"))
      print(model_type_use)
      print("ansamble7")
      time_taken <- system.time(my_ansamble(list7))
      print(time_taken)
      sink()
      print(time_taken)
  }
}