
#### Script for studying magnetic orientation of raptors, here REKI 
#### Script written by Filibert Heim, filibert.heim@posteo.de, after instructions from Steffen Oppel and hints from Will Scheider 

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 1. Preparations  --------
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# load required packages 
library(oce)
library(sf)
library(terra)
library(tidyverse)
select <- dplyr::select
filer <- dplyr::filter 
rename <- dplyr::rename

# read in data 
data <- read.csv(file = 'data/REKI_sample_locations10.csv') 

# check data structure 
head(data)
str(data)

# improve data structure 
data <- data %>% 
  mutate(timestamp = as.POSIXct(timestamp, format = "%Y-%m-%dT%H:%M:%SZ", tz = "UTC"), 
         migration = as.factor(if_else(MIGRATION == 'migrating', 'migratory', 'stationary')),
         long = as.numeric(sub("\\|.*", "", geometry)), 
         lat = as.numeric(sub(".*\\|", "", geometry))) %>% # also as.character
  select(-MIGRATION, -geometry)


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 2. Create a raster for the migration route of REKI's  --------
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

rast

######################################################
#### STEPS SUGGESTED BY STEFFEN ######################
######################################################

# 1. Ein 50 km raster über die bounding box des Zugweges erstellen, und dann für einen mittleren Tag während des Frühjahrszuges das magneticField für alle diese Rasterpunkte mit oce::magneticField die dec,inc, int Werte extrahieren
# 2. Für jeden Herbstzug (und Vogel), die inc, dec, und int Werte aller locations auf eine Toleranz runden welche der Wahrnehmungsfähigkeit entspricht (z.B. 0.5 für dec, die anderen müsste man nachschauen oder fragen).
# 3. Nach dem runden die unique() combinations aus inc, dec und int heraussuchen (geht zur Not mit paste, dann unique).
# 4. Das Raster über den Zugweg mit der gleichen sensitivität runden und auch die inc_dec, int Werte zusammen-pasten in einen Wert.
# 5. Für jede unique inc/dec/int combination alle locations aus dem raster rausfiltern die den gleichen (pasted) Wert von inc/dec/int haben. All diese Rasterzellen auf einer Karte einfärben.
# 6. Frühjahrszug des selben Vogels auf Karte plotten und checken ob die Route die eingefärbten Rasterzellen verlässt. Geht auch mit sf::st_within um festzustellen ob die locations innerhalb dieser zellen liegen.
# 7. Loop über alle Individuen machen und solange in die Berge gehen 



plot(x = data$long, y = data$lat)
