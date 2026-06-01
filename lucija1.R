rm(list=ls())

#Set path to the working directory containing:
#1. INTERMAGNET data set - to be taken from the INTERMAGNET web-site, and reformatted 
#accordingly - #reformatted original data enclosed
#2. TEC data set, as derived from the RINEX GPS observations using GPS TEC software - enclosed

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
library(AppliedPredictiveModeling)
library(DataExplorer)
library(corrplot)
library(glmnet)
library(caret)
library(car)
library(forecast)
library(relaimpo)
library(DALEX)
library(fitdistrplus)
library(monomvn)

oldw <- getOption("warn")
options(warn = -1)

#Data input: TEC and dTEC data derived from RINEX GPS observations taken at Darwin, NT, Australia
#using Gopi Seemala's GPS TEC software - days in 2014

iono1 <- matrix(nrow = 0, ncol = 3)

for(j in 1:12){ 
  month_file_name <- paste('darwin_data/data_all_', j, '.csv', sep = '')
  data_month <- as.data.frame(read.csv(month_file_name, header = TRUE, sep = ',', skip = 0))
  diono <- as.data.frame(cbind(data_month['total_seconds'], data_month['TEC'], data_month['TEC']))
  colnames(diono) <- c('time', 'TEC', 'dTEC')
  iono1 <- rbind(iono1, diono)
}

tec1 <- as.ts(iono1$TEC)

k <- as.data.frame(read.csv('Klobuchar.csv', header = TRUE, sep = ',', skip = 0))

cat(k[1])
