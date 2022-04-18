Travel Time: Raster-Based Approach
================
DIME - ieConnect
March 22, 2005

## Set up

``` r
## Load packages
library(rmarkdown)
library(dplyr)
```

    ## 
    ## Attaching package: 'dplyr'

    ## The following objects are masked from 'package:stats':
    ## 
    ##     filter, lag

    ## The following objects are masked from 'package:base':
    ## 
    ##     intersect, setdiff, setequal, union

``` r
library(sp)
library(sf)
```

    ## Linking to GEOS 3.9.1, GDAL 3.2.3, PROJ 7.2.1; sf_use_s2() is TRUE

``` r
library(raster)
```

    ## 
    ## Attaching package: 'raster'

    ## The following object is masked from 'package:dplyr':
    ## 
    ##     select

``` r
library(gdistance)
```

    ## Loading required package: igraph

    ## 
    ## Attaching package: 'igraph'

    ## The following object is masked from 'package:raster':
    ## 
    ##     union

    ## The following objects are masked from 'package:dplyr':
    ## 
    ##     as_data_frame, groups, union

    ## The following objects are masked from 'package:stats':
    ## 
    ##     decompose, spectrum

    ## The following object is masked from 'package:base':
    ## 
    ##     union

    ## Loading required package: Matrix

    ## 
    ## Attaching package: 'gdistance'

    ## The following object is masked from 'package:igraph':
    ## 
    ##     normalize

``` r
library(velox)
library(osmdata)
```

    ## Data (c) OpenStreetMap contributors, ODbL 1.0. https://www.openstreetmap.org/copyright

``` r
## Define projection
UTM_PROJ_EPSG <- 3968
```

## Create example points

Create example points to calculate travel time between

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
locs_sp <- spTransform(locs_sp, CRS(paste0("+init=epsg:", UTM_PROJ_EPSG)))
```
