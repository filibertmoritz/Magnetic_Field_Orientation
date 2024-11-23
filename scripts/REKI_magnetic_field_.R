
#### Script for studying magnetic orientation of raptors, here REKI 
#### Script written by Filibert Heim, filibert.heim@posteo.de, after instructions from Steffen Oppel and hints from Will Scheider 

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 1. Preparations  --------
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# load required packages 
library(oce)
library(sf)
library(terra)
# library(stars) this is another packages that might have advantages for combining raster and sf workflows
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

# create column which tells whether the current migration was the first autumn migration 
data <- data %>%
  group_by(id) %>%
  mutate(
    first_migration = cumsum(migration == 'migratory' & lag(migration, default = 'stationary') != 'migratory'),
    first_migration = ifelse(migration == 'migratory', first_migration, NA)) %>%
  ungroup() 

data <- data %>% 
  group_by(id) %>% 
  filter(migration == 'migratory') %>% 
  mutate(
    first_migration = if_else(year(timestamp) == year(min(timestamp, na.rm = TRUE)), 1, first_migration), 
    first_migration = if_else(year(timestamp) == year(min(timestamp, na.rm = TRUE)) + 1 & timestamp < as.POSIXct(paste0(year(min(timestamp, na.rm = TRUE)) + 1, '-07-15 00:00:00'), format = "%Y-%m-%d %H:%M:%S"), 2, first_migration)) %>%
  filter(first_migration %in% c(1,2)) %>% select(id, timestamp, first_migration) %>% 
  right_join(data %>% select(-first_migration), by = join_by(id, timestamp))

data <- data %>% mutate(first_migration = factor(case_when(first_migration == 1 ~ 'autumn', 
                                             first_migration == 2 ~ 'spring', 
                                            TRUE ~ as.character(first_migration))))

View(data)
# bring data in long format for easier plotting
# data <- data %>% pivot_longer(cols = c(inclination, declination, intensity),
#                              names_to = "Type",
#                              values_to = "Value")

summary(data$first_migration)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 2. Create a raster for the migration route of REKI's  --------
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# create binding box 
bbox <- terra::ext(range(data$long), range(data$lat))

# create raster with geographic CRS for Europe where grid cell size varies with location
rast <- rast(ext = bbox, resolution = 0.5, crs = "EPSG:4258")

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 3. Get magnetic field values for autumn migration  --------
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# calculate mean autumn migration date a
data <- data %>%
  group_by(id) %>%
  mutate(mean_autumn = if_else(first_migration == 'autumn', 
                               mean(timestamp[first_migration == 'autumn'], na.rm = TRUE), 
                               as.POSIXct(NA)))

# create df to extract magnetic field values for autumn migration for every bird 

birds <- unique(data$id) # get ids from all birds 
mean_autumn <- data %>% filter(first_migration == 'autumn') %>% group_by(id) %>% summarise(mean_autumn = mean(timestamp)) # get mean autumn migration dates 

df_rast <- as.data.frame(xyFromCell(rast, 1:ncell(rast))) # get the coordinated from raster
df_rast <- df_rast %>% rename(long = x, lat = y) # rename raster

df_rast <- data.frame(long = rep(df_rast$long, times = length(birds)), 
                      lat = rep(df_rast$lat, times = length(birds)), 
                      id = rep(birds, each = length(df_rast$long)))

df_rast <- df_rast %>% left_join(mean_autumn, by = join_by(id)) # add the mean_autumn migration dates 

for(i in 1:10){
  df_rast$inclination[i] <- magneticField(longitude = df_rast$long[i], latitude = df_rast$lat, time = df_rast$mean_autumn)$inclination[1]
  df_rast$declination[i] <- magneticField(longitude = df_rast$long[i], latitude = df_rast$lat, time = df_rast$mean_autumn)$declination[1]
  df_rast$intensity[i] <- magneticField(longitude = df_rast$long[i], latitude = df_rast$lat, time = df_rast$mean_autumn)$intensity[1]
}

df_rast <- df_rast %>%
  rowwise() %>%
  mutate(
    magnetic = list(magneticField(longitude = long, latitude = lat, time = mean_autumn)),
    inclination = magnetic$inclination,
    declination = magnetic$declination,
    intensity = magnetic$intensity) %>%
  select(-magnetic) %>%  # Remove the list column
  ungroup()

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



data %>% drop_na(first_migration) %>%
  ggplot()+ 
  geom_point(mapping = aes(x = timestamp, y = inclination, color = first_migration), size = 0.5) + 
  facet_wrap(~id, scales = 'free')
