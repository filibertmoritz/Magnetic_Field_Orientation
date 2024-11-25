
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
  select(-MIGRATION, -geometry, -inclination, -declination, -intensity)

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
# 3. Get magnetic field values for first spring migration in raster  --------
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# calculate mean autumn migration date a
data <- data %>%
  group_by(id) %>%
  mutate(mean_spring = if_else(first_migration == 'spring', 
                               mean(timestamp[first_migration == 'spring'], na.rm = TRUE), 
                               as.POSIXct(NA)))

# create df to extract magnetic field values for autumn migration for every bird 

birds <- unique(data$id) # get ids from all birds 
mean_spring <- data %>% filter(first_migration == 'spring') %>% group_by(id) %>% summarise(mean_spring = mean(timestamp)) # get mean autumn migration dates 

df_rast <- as.data.frame(xyFromCell(rast, 1:ncell(rast))) # get the coordinated from raster
df_rast <- df_rast %>% rename(long = x, lat = y) # rename raster

df_rast <- data.frame(long = rep(df_rast$long, times = length(birds)), 
                      lat = rep(df_rast$lat, times = length(birds)), 
                      id = rep(birds, each = nrow(df_rast)),
                      inclination = NA, 
                      declination = NA, 
                      intensity = NA)

df_rast <- df_rast %>% left_join(mean_spring, by = join_by(id)) # add the mean_spring migration dates 

df_rast$inclination <- magneticField(longitude = df_rast$long, latitude = df_rast$lat, time = df_rast$mean_spring)$inclination
df_rast$declination <- magneticField(longitude = df_rast$long, latitude = df_rast$lat, time = df_rast$mean_spring)$declination
df_rast$intensity <- magneticField(longitude = df_rast$long, latitude = df_rast$lat, time = df_rast$mean_spring)$intensity

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 3. Get magnetic field values for first autumn migration for all locations of each bird --------
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# filter for autumn migration
data_autumn <- data %>% 
  filter(first_migration == 'autumn') 

# get magnetic field values for all locations in first autumn migration
data_autumn$declination <- magneticField(longitude = data_autumn$long, latitude = data_autumn$lat, time = data_autumn$timestamp)$declination
data_autumn$inclination <- magneticField(longitude = data_autumn$long, latitude = data_autumn$lat, time = data_autumn$timestamp)$inclination
data_autumn$intensity <- magneticField(longitude = data_autumn$long, latitude = data_autumn$lat, time = data_autumn$timestamp)$intensity

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 4. Get information of grid cells where autumn migration values lied within spring migration  --------
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# round values which represents birds ability to sense magnetic field values, ATTENTION: only autumn migration magnetic field values are rounded, not spring migration
data_autumn$declination_max <- data_autumn$declination + 0.5 # declination sensitivity 0.5
data_autumn$declination_min <- data_autumn$declination - 0.5 
data_autumn$inclination_max <- data_autumn$inclination + 0.5 # inclination sensitivity 0.5
data_autumn$inclination_min <- data_autumn$inclination - 0.5 
data_autumn$intensity_max <- data_autumn$intensity + 200 # intensity sensitivity 200
data_autumn$intensity_min <- data_autumn$intensity - 200

# create one df of data_autumn and df_rast from spring migration for comparison
overlap <- df_rast %>% left_join(data_autumn %>% select(-long, -lat, -mean_spring), by = join_by(id), suffix = c('_spring', '_autumn')) %>% 
  select(-migration, -first_migration, -distance_to_next, -speed) %>% 
  mutate(match_inclination = inclination_spring >= inclination_min & inclination_spring <= inclination_max, 
         match_declination = declination_spring >= declination_min & declination_spring <= declination_max,
         match_intensity = intensity_spring >= intensity_min & intensity_spring <= intensity_max, 
         match_all = match_inclination & match_declination & match_intensity) %>% 
  select(-match_inclination, -match_declination, -match_intensity) %>% 
  filter(match_all == TRUE) 

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 5. Plot data --------
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# get information on raster 
str(rast)
dim(rast)

# create a plot 
overlap %>% 
  filter(id == 'REKI_104') %>% 
  ggplot() + 
  geom_point(mapping = aes(x = long, y = lat), size = 2) +
  coord_equal()


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
