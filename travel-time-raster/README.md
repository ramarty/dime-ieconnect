# Computing Travel Time using Raster Based Approach

Functions to facilitate estimating travel time, shortest paths, and market access using a raster-based approach. Here, a polyline of a road network is rasterized and converted into a transition object, enabling travel time and shortest paths to be calculated across the area.

## Example
See example [here](https://github.com/ramarty/dime-ieconnect/blob/main/travel-time-raster/example/travel-time-raster-example.md).

## Load functions
```r
source("https://raw.githubusercontent.com/ramarty/dime-ieconnect/main/travel-time-raster/travel-time-raster.R")
```

## Functions
```r
#### Rasterize roads
# Rasterize a road network. Takes the network and, for each pixel within the raster, assigns the fastest speed. 
rasterize_roads(road_sdf, 
                speed_var, 
                extent_sdf = NULL,
                pixel_size_km = 0.5,
                walking_speed = 5,
                pixel_var = "time_to_cross",
                restrict_to_extent = F)

#### Make travel time matrix
# Calculates the travel time for every point in `points_sdf` to all other points in `points_sdf`. Returns a matrix with the following variables: 
# --orig_[uid_name]: uid of origin location
# --dest_[uid_name]: uid of destination location
# --travel_time: travel time between origin-destination (in hours)
# --distance_meters: distance (in meters) between origin-destination
make_tt_matrix(points_sdf,
               uid_name,
               cost_t)


#### Calculate market access
# Relying on a travel time matrix computed using `make_tt_matrix` and a measure of market size (e.g., population), computes market access for each location
calc_ma(tt_df,
        orig_uid_var,
        market_var,
        travel_cost_var,
        theta,
        exclude_km)
```




