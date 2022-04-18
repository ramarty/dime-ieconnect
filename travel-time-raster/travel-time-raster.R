# Functions to faciliate estimating travel times

library(dplyr)
library(sp)
library(raster)
library(gdistance)
library(velox)
library(rgeos)
library(stringr)

rasterize_roads <- function(road_sdf, 
                            speed_var, 
                            extent_sdf = NULL,
                            pixel_size_km = 0.5,
                            walking_speed = 5,
                            pixel_var = "time_to_cross",
                            restrict_to_extent = F){
  # DESCRIPTION: Rasterize a road network. Takes the network and, for each pixel
  # within the raster, assigns the fastest speed. 
  # ARGS:
  #   road_sdf: Spatial dataframe of roads
  #   speed_var: (character) name of variable in spatial dataframe that corresponds to speeds 
  #   extent: Spatial object where extent should be derived for creating raster. The function
  #           will create a raster using this bounding box. If NULL, uses extent of road_sdf
  #   walking_speed: Walking speed (speed assigned to areas not roads). 
  #   pixel_var: What the pixel value should represent. 
  #              (a) time_to_cross: gives the time to cross each pixel
  #              (b) speed: speed of pixel
  #   restrict_to_extent: If TRUE, areas outside of extent_df are set to 0 (so 
  #                       can't travel on the area). E.g., if extent is a country,
  #                       can't travel outside of the country
  # NOTE: Assumes both road_sdf and extent (1) have the same CRS and (2) are projected
  #       (so map units are meters)
  
  #### Prep Road Obj 
  
  ## Prep speed var
  road_sdf$speed_var_temp <- road_sdf[[speed_var]]
  road_sdf <- road_sdf[!is.na(road_sdf$speed_var_temp),]
  road_sdf <- road_sdf[!(road_sdf$speed_var_temp %in% 0),]
  
  ## Sort by Speed
  # If multiple polylines interesect with a cell, velox uses the last polygon from
  # the spatial polygons dataframe. Consequently, we sort by speeds from slowest to
  # fastest so that velox uses the fastest speed.
  road_sdf <- road_sdf[order(road_sdf$speed_var_temp, decreasing=F),] 
  
  #### Prep extent object
  if(is.null(extent_sdf)) extent_sdf <- road_sdf
  
  #### Make blank raster
  r <- raster(xmn=extent_sdf@bbox[1,1], 
              xmx=extent_sdf@bbox[1,2], 
              ymn=extent_sdf@bbox[2,1], 
              ymx=extent_sdf@bbox[2,2], 
              crs=crs(road_sdf), 
              resolution = pixel_size_km*1000)
  
  r <- raster(ext = extent_sdf %>% extent(), 
              crs=crs(road_sdf), 
              resolution = pixel_size_km*1000)
  
  #### Rasterize Speeds
  roads_r <- r
  roads_r[] <- 0
  roads_r_vx <- velox(roads_r)
  roads_r_vx$rasterize(road_sdf, field="speed_var_temp", background=walking_speed) # background should be walking speed (5km/hr); https://en.wikipedia.org/wiki/Preferred_walking_speed
  roads_r <- roads_r_vx$as.RasterLayer()
  
  #### Convert from speed to time to cross
  if(pixel_var %in% "time_to_cross") roads_r[] <- pixel_size_km / roads_r[]
  
  #### Mask
  if(restrict_to_extent %in% T){
    
    roads_r <- raster::mask(roads_r, extent_sdf)
    
    roads_r[][is.na(roads_r[])] <- (pixel_size_km/walking_speed) * 100
    
  }
  
  return(roads_r)
}

make_tt_matrix <- function(points_sdf,
                           uid_name,
                           cost_t){
  
  # DESCRIPTION: Calculates the travel time for every point in `points_sdf` to
  # all other points in `points_sdf`. Returns a matrix with the following variables
  # orig_[uid_name]: uid of origin location
  # dest_[uid_name]: uid of destination location
  # travel_time: travel time between origin-destination (in hours)
  # distance_meters: distance (in meters) between origin-destination
  
  ## Give uid variable
  points_sdf$uid <- points_sdf[[uid_name]]
  
  tt_df <- lapply(1:nrow(points_sdf), function(i){
    if((i %% 10) == 0) print(i)
    
    tt <- costDistance(cost_t,
                       points_sdf[i,],
                       points_sdf) %>% as.numeric()
    
    df_out <- data.frame(dest_uid    = points_sdf$uid,
                         travel_time = tt)
    df_out$distance_meters <- gDistance(points_sdf[i,], points_sdf, byid = T) %>% as.numeric()
    
    df_out$orig_uid <- points_sdf$uid[i]
    
    ## Replace "uid" with [uid_name]
    names(df_out) <- names(df_out) %>% str_replace_all("uid", uid_name)
    return(df_out)
  }) %>% bind_rows()
  
  return(tt_df)
}

calc_ma <- function(tt_df,
                    orig_uid_var,
                    market_var,
                    travel_cost_var,
                    theta,
                    exclude_km){
  # DESCRIPTION: Using a travel cost matrix that includes market size (eg, population),
  # computes market access.
  # ARGS:
  #  tt_df: Travel time matrix dataframe
  #  orig_uid_var: Name of origin uid variable
  #  market_var: Name of variable to use as market size
  #  travel_cost_var: Name of variable to use as travel cost
  #  theta: Cost parameter (theta)
  #  exclude_km: When calculating market access, exclude units within this distance
  
  ## Prep Variables
  tt_df$orig_uid    <- tt_df[[orig_uid_var]]
  tt_df$market_var  <- tt_df[[market_var]]
  tt_df$travel_cost <- tt_df[[travel_cost_var]]
  
  ## Compute Market Access
  MA_df <- tt_df[tt_df$distance_meters > exclude_km*1000,] %>% 
    group_by(orig_uid) %>%
    dplyr::summarise(MA_VAR  = sum(market_var / (travel_cost^theta)))
  
  ## Names for variable name
  travel_cost_name <- ""
  if(travel_cost_var %in% "travel_time") travel_cost_name <- "tt"
  if(grepl("distance", travel_cost_var)) travel_cost_name <- "dist"
  theta_name <- theta %>% as.character %>% str_replace_all("\\.", "_")
  if(exclude_km %in% 0) exclude_name <- ""
  if(exclude_km > 0) exclude_name <- paste0("_exclude", exclude_km, "km")
  
  names(MA_df)[names(MA_df) %in% "MA_VAR"] <- paste0("MA_",
                                                     travel_cost_name,
                                                     "_",
                                                     market_var,
                                                     "_theta",
                                                     theta_name,
                                                     exclude_name)
  
  return(MA_df)
  
}