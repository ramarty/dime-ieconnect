# Create example lines

library(dplyr)
library(sp)
library(sf)
library(raster)
library(gdistance)
library(velox)
library(osmdata)

UTM_PROJ_EPSG <- 3968

# Example locations for travel time --------------------------------------------
locs_sp <- bind_rows(
  data.frame(name = "national cathedral",
             lat = 38.930589,
             lon = -77.070803),
  
  data.frame(name = "capitol",
             lat = 38.889867,
             lon = -77.009064),
  
  data.frame(name = "world bank",
             lat = 38.898987,
             lon = -77.042402),
  
  data.frame(name = "takoma metro",
             lat = 38.975616,
             lon = -77.017906)
)

## Add unique ID
locs_sp$id <- 1:nrow(locs_sp)

## Spatially define
coordinates(locs_sp) <- ~lon+lat
crs(locs_sp) <- CRS("+init=epsg:4326")

## Project
locs_sp <- spTransform(locs_sp, CRS(paste0("+init=epsg:", UTM_PROJ_EPSG)))

# Grab and clean road data from OSM --------------------------------------------
#### Import OSM road data
q <- opq(bbox = 'washington dc') %>%
  add_osm_feature(key = 'highway', value = c('trunk',
                                             'primary')) %>%
  osmdata_sf()

roads_sf <- q$osm_lines

#### Assign speed limits
roads_sf$speed_kmhr <- NA
roads_sf$speed_kmhr[roads_sf$highway %in% "primary"] <- 30
roads_sf$speed_kmhr[roads_sf$highway %in% "trunk"]   <- 60

roads_sf <- st_transform(roads_sf, 3968)

#### To SP
roads_sp <- roads_sf %>% as("Spatial")

# Make transition object -------------------------------------------------------
#### Rasterize roads
roads_r <- rasterize_roads(road_sdf = roads_sp,
                           speed_var = "speed_kmhr",
                           pixel_size_km = 0.1)

#### Make transition object
cost_t <- transition(roads_r, function(x) 1/mean(x), directions=8)

# Travel times+routes from one location to all others --------------------------
wb_sp       <- locs_sp[locs_sp$name == "world bank",]
otherloc_sp <- locs_sp[locs_sp$name != "world bank",]

spath_sp <- shortestPath(cost_t, 
                         wb_sp, 
                         otherloc_sp, 
                         output="SpatialLines")
spath_sp$travel_time_hr <- costDistance(cost_t,
                                        wb_sp,
                                        otherloc_sp) %>% as.numeric()
spath_sp$name <- otherloc_sp$name

#### Check data
plot(roads_r)
plot(spath_sp, add = T)
head(spath_sp)

# Travel time matrix -----------------------------------------------------------
make_tt_matrix(points_sdf = locs_sp,
               uid_name   = "id",
               cost_t     = cost_t)











