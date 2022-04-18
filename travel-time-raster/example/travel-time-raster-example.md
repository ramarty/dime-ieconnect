Compute Travel Time: Raster-Based Approach
================

-   [Prep Points and and Roads Data](#prep-points-and-and-roads-data)
    -   [Set up](#set-up)
    -   [Create example points](#create-example-points)
-   [Calculate travel time, shortest path and market
    access](#calculate-travel-time-shortest-path-and-market-access)
    -   [Rasterize roads and make transition
        object](#rasterize-roads-and-make-transition-object)
    -   [Travel times & shortest paths from one location to all
        others](#travel-times--shortest-paths-from-one-location-to-all-others)
    -   [Travel time matrix](#travel-time-matrix)
    -   [Market access](#market-access)

This document describes computing travel time using a raster-based
approach. “Raster-based” means we convert poylines into a raster, where
each pixel in the raster indicates the time it takes to travel across
the raster (for example, for a 1km pixel, if a 60km/h road passes
through the pixel, it would take 1 minute to pass through the pixel). We
can then convert this raster into a “transition object”, which allows
easily computing travel times and shortest paths between any two
locations relying on the [gdistance
package](https://cran.microsoft.com/snapshot/2014-12-09/web/packages/gdistance/vignettes/gdistance.pdf).

The first part of the document creates and prepares (1) point data for
locations in Washington, DC where we want to calculate travel times
between and (2) road data for Washington, DC.

The second part of the document creates the transition object, then uses
that object to compute travel times, shortest paths, and a measure of
market access.

Instead of a raster-based approach, one could convert the road polyline
into graph, and compute travel times between nodes in the graph. The
benefit of a raster-based approach is that it doesn’t require converting
the road network polyline into a graph, which in some instances can
require a fair bit of work. However, see
[here](https://r-spatial.org/r/2019/09/26/spatial-networks.html) for an
example of a using this graph/network-based approach.

# Prep Points and and Roads Data

Here we: 1. Prep a `SpatialPoints` object of locations to calculate
travel times between 2. Prep a `SpatialPolyline` object of roads

## Set up

Here, we (1) load packages, (2) import travel time functions, and (3)
define the projection to use on the data (the travel time functions
assume the data are projected).

``` r
## Load packages
library(rmarkdown)
library(dplyr)
library(sp)
library(sf)
library(rgeos)
library(raster)
library(gdistance)
library(velox)
library(osmdata)
library(stringr)
library(leaflet)

## Travel Time & MA Functions
source("https://raw.githubusercontent.com/ramarty/dime-ieconnect/main/travel-time-raster/travel-time-raster.R")

## Define projection
UTM_PROJ <- "+init=epsg:3968"
```

## Create example points

Create example points to calculate travel time between.

``` r
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
locs_sp <- spTransform(locs_sp, CRS(UTM_PROJ))
```

``` r
leaflet() %>%
  addTiles() %>%
  addCircles(data = locs_sp %>% spTransform("+init=epsg:4326"), 
             popup = ~name,
             weight = 15,
             opacity = 1)
```

![](travel-time-raster-example_files/figure-gfm/unnamed-chunk-3-1.png)<!-- -->
\## Grab and clean road data from OSM Here, we grab road network data
from OpenStreetMaps for Washington, DC. The code for calculating travel
times assumes the spatial data is projected.

``` r
#### Import OSM road data
q <- opq(bbox = 'washington dc') %>%
  add_osm_feature(key = 'highway', value = c('trunk',
                                             'primary')) %>%
  osmdata_sp()

roads_sp <- q$osm_lines

#### Assign speed limits
roads_sp$speed_kmhr <- NA
roads_sp$speed_kmhr[roads_sp$highway %in% "primary"] <- 30
roads_sp$speed_kmhr[roads_sp$highway %in% "trunk"]   <- 60

## Project
roads_sp <- spTransform(roads_sp, UTM_PROJ)
```

# Calculate travel time, shortest path and market access

The below code illustrates using the road network data to calculate
travel time between different locations. We: 1. Rasterize roads and make
a transition object. The transition object allows us to pick any two
points on the surface and easily get the travel time and shortest path
2. Compute travel time and shortest paths from one location to all
locations 3. Create a travel time matrix that includes the travel time
between all locations 4. Compute a measure of market access for each
location

## Rasterize roads and make transition object

``` r
#### Rasterize roads
roads_r <- rasterize_roads(road_sdf = roads_sp,
                           speed_var = "speed_kmhr",
                           pixel_size_km = 0.1)

#### Make transition object
cost_t <- transition(roads_r, function(x) 1/mean(x), directions=8)
```

``` r
plot(roads_r)
```

![](travel-time-raster-example_files/figure-gfm/unnamed-chunk-6-1.png)<!-- -->

## Travel times & shortest paths from one location to all others

We calculate the travel time from the World Bank to the other locations
and create polylines of the shortest paths.

``` r
## Make separate spatial dataframes for (1) origin and (2) destinations 
o_sp  <- locs_sp[locs_sp$name == "world bank",]
ds_sp <- locs_sp[locs_sp$name != "world bank",]

## Create polylines of shortest paths
od_paths_sp <- shortestPath(cost_t, 
                         o_sp, 
                         ds_sp, 
                         output="SpatialLines")

## Add Travel time
od_paths_sp$travel_time_hr <- costDistance(cost_t,
                                        o_sp,
                                        ds_sp) %>% as.numeric()

## Add o-d names to dataframe
od_paths_sp$o_name <- "world bank"
od_paths_sp$d_name <- ds_sp$d_name
```

``` r
head(od_paths_sp)
```

``` r
plot(roads_r)
plot(od_paths_sp, add = T)
```

![](travel-time-raster-example_files/figure-gfm/unnamed-chunk-9-1.png)<!-- -->

## Travel time matrix

Here, we create a travel time matrix that computes the travel time from
each location to all other locations.

``` r
tt_df <- make_tt_matrix(points_sdf = locs_sp,
                        uid_name   = "id",
                        cost_t     = cost_t)
```

``` r
print(tt_df)
```

## Market access

Here, we compute a measure of market access. We rely on the travel time
matrix computed above. For the market access measure, we use population
(for example, this could be the population of the neighborhood in
Washington, DC; we just randomly generate this for illustration
purposes). Market access is computed as (where
![o](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;o "o")
refers to origin and
![d](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;d "d")
destination):

![MA\_{o} = \\sum\_{d,\~d \\neq o} pop_d \\times tt\_{o,d}^{-\\theta}](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;MA_%7Bo%7D%20%3D%20%5Csum_%7Bd%2C~d%20%5Cneq%20o%7D%20pop_d%20%5Ctimes%20tt_%7Bo%2Cd%7D%5E%7B-%5Ctheta%7D "MA_{o} = \sum_{d,~d \neq o} pop_d \times tt_{o,d}^{-\theta}")

``` r
set.seed(42)
locs_sp$population <- runif(nrow(locs_sp))*10000

tt_data_df <- tt_df %>%
  left_join(locs_sp@data, by = c("dest_id" = "id"))

ma_df <- calc_ma(tt_df = tt_data_df,
                 orig_uid_var = "orig_id",
                 market_var = "population",
                 travel_cost_var = "travel_time",
                 theta = 3,
                 exclude_km = 0)
```

``` r
print(ma_df)
```
