rm(list=ls())

# Declaration of the R libraries utilised

library(ggplot2)
library(sf)
library(ggspatial)
library(tidyverse)
library(prettymapr)

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

stations <- data.frame(
  place=c('Darwin, NT', 'Kakadu, NT'),
  lon=c(131.1327361,132.47),
  lat=c(-12.8437111,-12.69),
  nudge_y=c(-0.1,0.1),
  nudge_x=c(-0.2,0.2) 
)

stations_sf <- st_as_sf(stations,coords=c('lon','lat'),crs=4326)

map_plot<-ggplot() +
  annotation_map_tile(type='osm',zoomin=0) +
  geom_sf(data=stations_sf,color='red',size=3) +
  geom_sf_text(
    data=stations_sf, 
    aes(label=place), 
    nudge_x=stations$nudge_x,
    nudge_y=stations$nudge_y,
    fontface='bold',color='black',
    size=4
  ) +
  coord_sf(
    xlim=c(130.5,133.0), 
    ylim=c(-13.5,-12.0),
    crs=4326,
    datum=4326
  ) +
  theme_minimal() +
  labs(
    title='Map of observing stations',
    x='Longitude',
    y='Latitude'
  )

dir.create(file.path(paste('TEC_',limited,sep='')),showWarnings=FALSE)
ggsave(
  filename=paste('TEC_',limited,'/Observing_stations_map.pdf_',limited,'.pdf',sep=''), 
  plot=map_plot,
  width=8,         # Width of the PDF
  height=6,        # Height of the PDF
  units='in',      # Units (inches)
  dpi=600          # Ensures the OSM background tiles are high-resolution
)