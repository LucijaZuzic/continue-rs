rm(list=ls())

# Declaration of the R libraries utilised

library(tidyverse)

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

# Process all days of the year (1 to 365)

results<-lapply(1:365,function(DOY) {
  
  # Construct the file path
  
  file_path<-paste0('GPS navigation messages 2014/brdc',DOY,'0.14n')
  
  # Check if file exists
  
  if (!file.exists(file_path)) return(NULL)
  
  # Read the first 5 lines of the file
  
  head<-readLines(file_path,n=5,warn=FALSE)
  
  # Extract lines 4 and 5
  
  alphas_str<-head[4]
  betas_str<-head[5]
  
  # Clean, split, replace 'D' with 'E', and convert to numeric
  # trimws() removes leading/trailing spaces
  # strsplit(..., '\\s+') splits by any amount of whitespace
  # gsub() replaces the Fortran 'D' scientific notation with standard 'E'
  
  alphas<-as.numeric(gsub('D','E',strsplit(trimws(alphas_str),'\\s+')[[1]]))
  betas<-as.numeric(gsub('D','E',strsplit(trimws(betas_str),'\\s+')[[1]]))
  
  # Return as a 1-row data frame
  
  data.frame(
    DOY=DOY,
    A1=alphas[1],A2=alphas[2],A3=alphas[3],A4=alphas[4],
    B1=betas[1],B2=betas[2],B3=betas[3],B4=betas[4]
  )
})

# Combine the list of 1-row data frames into a single data frame

new_df<-do.call(rbind,results)

# Write to CSV (row.names=FALSE to remove the index column)

write.csv(new_df,'Klobuchar.csv',row.names=FALSE)