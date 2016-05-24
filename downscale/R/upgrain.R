################################################################################
# 
# upgrain.R
# Version 1.4
# 25/02/2016
#
# Updates:
#   25/02/2016: Bug on extent in standardised rasters fixed
#               Added option to save all rasters
#               Change the method of assigning raster values to 'setValues'
#   05/10/2015: Bug on number of scales fixed
#               'All_Presences' changed to 'All_occurrences'
#               Plotting now optional
#   30/04/2015: Threshold selection added
#   14/04/2015: plotting functions added
#               help file updated
#
# Upgrains atlas data to larger grain sizes. All raster files are then
# standardised to have the same extent as the extent of the maximum grain size.
# All new cells are given a value of 0.
#
# data = either a raster file of presence-absence atlas data, or a data frame of 
#        sampled cells with longitude, latitude and presence-absence.
#
# cell.width = if data is a data frame, the cell widths of sampled cells.
#
# scales = the number of cells to upgrain. Upgraining will happen by factors of
#          2 - ie if scales = 3, the atlas data will be aggregated in 2x2 cells,
#          4x4 cells and 8x8 cells.
#
################################################################################

upgrain <- function(atlas.data,
                    cell.width = NULL,
                    scales,
                    threshold = NULL,
                    method = "Gain_Equals_Loss",
                    plot = TRUE,
                    return.rasters = FALSE) {
  
  ### Error checking: either a threshold is given or a method selected, not both
  if((is.null(threshold) == FALSE) & (is.null(method) == FALSE)) {
    stop("Both a threshold and a model have been given. Please set one to NULL")
  }
  
  if(is.null(method)) {
    ### Error checking: threshold is between 0 and 1
    if(min(threshold) < 0){
      stop("Threshold value must be between 0 and 1")
    }
    if(max(threshold) > 1){
      stop("Threshold value must be between 0 and 1")
    }
  }
  
  if(is.null(threshold)) {
    ### Error checking: only one method is selected
    if(sum(method == c("Sampled_Only",
                       "All_Sampled",
                       "All_Occurrences",
                       "Gain_Equals_Loss")) > 1){
      stop("Method is not one of the options")
    }
    ### Error checking: method is one of the selections
    if(sum(method == c("Sampled_Only",
                       "All_Sampled",
                       "All_Occurrences",
                       "Gain_Equals_Loss")) == 0){
      stop("Method is not one of the options")
    }
  }
  
  ### Error checking: scales
  if(scales < 2){
    stop("Scales must be >1 (at least three grain sizes needed for downscaling")
  }
  
  ### Error checking: if data frame needs cell width
  if(is.data.frame(atlas.data) == "TRUE") {
    if(is.null(cell.width)) {
      stop("If data is data.frame cell.width is required")
    }
  }
  
  ################################################################################
  ### data storage
  original <- data.frame("Cell.area" = rep(NA, scales + 1), 
                         "Extent"  = rep(NA, scales + 1),
                         "Occupancy" = rep(NA, scales + 1))
  extended <- data.frame("Cell.area" = rep(NA, scales + 1), 
                         "Extent"  = rep(NA, scales + 1),
                         "Occupancy" = rep(NA, scales + 1))
  
  ### data manipulation
  if(is.data.frame(atlas.data) == "TRUE") {
    longitude <- atlas.data[, 1]
    latitude <- atlas.data[, 2]
    presence <- atlas.data[, 3]
    presence[presence > 0] <- 1
    shapefile <- sp::SpatialPointsDataFrame(coords = data.frame(lon = longitude,
                                                                lat = latitude),
                                            data = data.frame(presence = presence)) 
    atlas_raster <- raster::raster(ymn = min(latitude) - (cell.width / 2),
                                   ymx = max(latitude) + (cell.width / 2),
                                   xmn = min(longitude) - (cell.width / 2), 
                                   xmx = max(longitude) + (cell.width / 2),
                                   resolution = cell.width)
    atlas_raster <- raster::rasterize(shapefile, atlas_raster)
    atlas_raster <- raster::dropLayer(atlas_raster, 1)
  }
  
  if(class(atlas.data) == "RasterLayer") {
    atlas_raster <- atlas.data@data@values
    atlas_raster[atlas_raster > 0] <- 1
    atlas_raster <- setValues(atlas.data, atlas_raster)
  }
  
  ### largest scale (all other layers are extended to equal this raster)
  max_raster <- raster::aggregate(atlas_raster, (2 ^ scales), fun = max)
  atlas_raster_extend <- raster::extend(atlas_raster, 
                                        raster::extent(max_raster))
  
  ###############################################################
  ########## Set new atlas data to selected threshold (atlas_thresh)
  
  ### create boundary raster
  boundary_raster <- atlas_raster
  boundary_raster@data@values[!is.na(boundary_raster@data@values)] <- 1
  boundary_raster <- raster::aggregate(boundary_raster, (2 ^ scales), fun = sum)  
  boundary_raster@data@values <- boundary_raster@data@values / ((2 ^ scales) ^ 2)
  
  ### run for threshold, "Sampled_Only" and "All_Sampled"
  if(is.null(threshold)) {
    if(method == "Sampled_Only") {
      threshold <- 1
    }
    if(method == "All_Sampled") {
      threshold <- 0
    }
  }
  if(is.null(threshold) == FALSE) {
    ### set extent of maximum grain size raster to threshold (max_raster_thresh)
    max_raster_thresh <- max_raster
    max_raster_thresh@data@values[boundary_raster@data@values < threshold] <- NA
    max_raster_thresh@data@values[max_raster_thresh@data@values == 0] <- 1
    max_raster_thresh <- raster::disaggregate(max_raster_thresh, 
                                              (2 ^ scales), fun = max)
    
    ### extend and crop atlas raster to new extent (atlas_thresh)
    atlas_thresh <- atlas_raster
    atlas_thresh <- raster::extend(atlas_thresh, 
                                   raster::extent(max_raster_thresh))
    atlas_thresh@data@values[atlas_thresh@data@values == 1] <- 2
    atlas_thresh@data@values[atlas_thresh@data@values == 0] <- 1
    atlas_thresh@data@values[is.na(atlas_thresh@data@values)] <- 0
    atlas_thresh <- raster::overlay(max_raster_thresh, atlas_thresh, fun = sum)
    atlas_thresh@data@values[atlas_thresh@data@values == 1] <- 0
    atlas_thresh@data@values[atlas_thresh@data@values == 2] <- 0
    atlas_thresh@data@values[atlas_thresh@data@values == 3] <- 1
  }
  
  ### run for "All_Occurrences" or "Gain_Equals_Loss")
  if(is.null(threshold)) {
    thresholds <- seq(0, 1, 0.01)
    
    ### Loop through thresholds
    land <- data.frame(Threshold = thresholds, 
                       SampledExcluded = rep(NA, length(thresholds)),
                       SampledIncluded = NA,
                       UnsampledAdded = NA, 
                       Extent = NA,
                       OccurrencesExcluded = NA)
    
    atlas_boundary <- atlas_raster_extend
    atlas_boundary@data@values[atlas_boundary@data@values == 0] <- 1
    atlas_boundary@data@values[is.na(atlas_boundary@data@values)] <- 0
    
    for(j in 1:length(thresholds)){
      atlas_boundary_thresh <- max_raster
      atlas_boundary_thresh@data@values[boundary_raster@data@values <
                                          thresholds[j]] <- NA
      atlas_boundary_thresh@data@values[atlas_boundary_thresh@data@values >=
                                          0] <- 2
      atlas_boundary_thresh@data@values[is.na(
        atlas_boundary_thresh@data@values)] <- 0
      atlas_boundary_thresh <- raster::disaggregate(atlas_boundary_thresh,
                                                    (2 ^ scales),
                                                    fun = max)
      both <- raster::overlay(atlas_boundary_thresh, atlas_boundary, fun = sum)
      land[j, "SampledExcluded"] <- sum(both@data@values == 1, na.rm = TRUE)
      land[j, "UnsampledAdded"] <- sum(both@data@values == 2, na.rm = TRUE)
      land[j, "SampledIncluded"] <- sum(both@data@values == 3, na.rm = TRUE)
      land[j, "Extent"] <- sum(atlas_boundary_thresh@data@values == 2, 
                               na.rm = TRUE)
      
      both <- raster::overlay(atlas_boundary_thresh, atlas_raster_extend,
                              fun = sum)
      land[j, "OccurrencesExcluded"] <- round(sum(both@data@values == 1, 
                                                  na.rm = TRUE) /
                                                sum(atlas_raster@data@values == 1,
                                                    na.rm = TRUE), 3)
    }
    Gain_loss.thresh <- thresholds[which.min(abs(land[, "Extent"] - 
                                                   sum(atlas_boundary@data@values, na.rm = TRUE)))]
    Presence.thresh <- thresholds[max(which(land[, "OccurrencesExcluded"]== 0))]
  }
  
  if(is.null(threshold)) {    
    if(method == "All_Occurrences") {
      threshold <- Presence.thresh
    }
    if(method == "Gain_Equals_Loss") {
      threshold <- Gain_loss.thresh
    }
    ### set extent of maximum grain size raster to threshold (max_raster_thresh)
    max_raster_thresh <- max_raster
    max_raster_thresh@data@values[boundary_raster@data@values < threshold] <- NA
    max_raster_thresh@data@values[max_raster_thresh@data@values == 0] <- 1
    max_raster_thresh <- raster::disaggregate(max_raster_thresh,
                                              (2 ^ scales), fun = max)
    
    ### extend and crop atlas raster to new extent (atlas_thresh)
    atlas_thresh <- atlas_raster
    atlas_thresh <- raster::extend(atlas_thresh, 
                                   raster::extent(max_raster))
    atlas_thresh@data@values[atlas_thresh@data@values == 1] <- 2
    atlas_thresh@data@values[atlas_thresh@data@values == 0] <- 1
    atlas_thresh@data@values[is.na(atlas_thresh@data@values)] <- 0
    atlas_thresh <- raster::overlay(max_raster_thresh, atlas_thresh, fun = sum)
    atlas_thresh@data@values[atlas_thresh@data@values == 1] <- 0
    atlas_thresh@data@values[atlas_thresh@data@values == 2] <- 0
    atlas_thresh@data@values[atlas_thresh@data@values == 3] <- 1
  }
  
  ####################################################################
  #### calculate occupancy at all grain sizes in original atlas data
  original[1, "Cell.area"] <- raster::res(atlas_raster)[1] ^ 2
  original[1, "Extent"] <- sum(!is.na(atlas_raster@data@values)) * 
    original[1, "Cell.area"]
  original[1, "Occupancy"] <- sum(atlas_raster@data@values == 1, na.rm = TRUE) /
    sum(!is.na(atlas_raster@data@values))
  
  for(i in 1:scales) {
    scaled_raster <- raster::aggregate(atlas_raster, 2 ^ i, fun = max)
    original[i + 1, "Cell.area"] <- raster::res(scaled_raster)[1] ^ 2
    original[i + 1, "Extent"] <- sum(!is.na(scaled_raster@data@values)) * 
      original[i + 1, "Cell.area"]
    original[i + 1, "Occupancy"] <- sum(scaled_raster@data@values == 1,
                                        na.rm = TRUE) /
      sum(!is.na(scaled_raster@data@values))
  }
  
  ####################################################################
  #### calculate occupancy at all grain sizes in standardised atlas data
  extended[1, "Cell.area"] <- raster::res(atlas_thresh)[1] ^ 2
  extended[1, "Extent"] <- sum(!is.na(atlas_thresh@data@values)) * 
    extended[1, "Cell.area"]
  extended[1, "Occupancy"] <- sum(atlas_thresh@data@values == 1, na.rm = TRUE) /
    sum(!is.na(atlas_thresh@data@values))
  
  for(i in 1:scales) {
    scaled_raster <- raster::aggregate(atlas_thresh, 2 ^ i, fun = max)
    if(return.rasters == TRUE) {
      assign(paste("scaled_raster", i, sep = "_"), scaled_raster)
    }
    extended[i + 1, "Cell.area"] <- raster::res(scaled_raster)[1] ^ 2
    extended[i + 1, "Extent"] <- sum(!is.na(scaled_raster@data@values)) *
      extended[i + 1, "Cell.area"]
    extended[i + 1, "Occupancy"] <- sum(scaled_raster@data@values == 1,
                                        na.rm = TRUE) /
      sum(!is.na(scaled_raster@data@values))
  }
  extended[, "Cell.area"] <- original[, "Cell.area"]
  
  ### plotting
  if(plot == TRUE) {  
    ### set par values for plotting
    par.original <- par()
    par.original <- list(mfrow = par.original$mfrow, mar = par.original$mar)
    par(mfrow = c(2, ceiling((scales + 2) / 2)), mar = c(3, 3, 3.5, 3))
    
    ########### Plot 1 - original atlas data
    plot(atlas_raster_extend,
         col = c("white", "white"),
         legend = FALSE, 
         axes = FALSE,
         main = paste("Original atlas data:\n cell area = ", 
                      original[1, "Cell.area"], sep = ""))
    plot(atlas_raster,
         col = c("white", "red"),
         legend = FALSE, 
         axes = FALSE,
         colNA = "dark grey",
         add = TRUE)
    
    ########### Plot 2 - standardised atlas data - atlas grain size
    plot(atlas_thresh,
         col = c("white", "red"),
         legend = FALSE, 
         axes = FALSE,
         colNA = "dark grey",
         main = paste("Standardised atlas data:\n cell area = ", 
                      original[1, "Cell.area"], sep = ""))
    
    for(i in 1:scales) {
      scaled_raster <- raster::aggregate(atlas_thresh, 2 ^ i, fun = max)
      
      ########### Plots 3 ... - standardised atlas data - all other grain sizes
      plot(scaled_raster,
           col = c("white", "red"),
           colNA = "dark grey",
           legend = FALSE, 
           axes = FALSE,
           main = paste("Standardised atlas data:\n cell area = ", 
                        original[i + 1, "Cell.area"], sep = ""))    
    }
    
    ### revert par to original values
    par(mfrow = par.original$mfrow, mar = par.original$mar)
  }
  
  output <- list(threshold = threshold,
                 extent.stand = extended[1, "Extent"],
                 occupancy.stand = extended,
                 occupancy.orig = original,
                 atlas.raster.stand = atlas_thresh)
  if(return.rasters == TRUE) {
    scaled.rasters <- list()
    for(i in 1:scales) {
      scaled.rasters[i] <- get(paste("scaled_raster", i, sep = "_"))
    }
    names(scaled.rasters) <- paste("scaled.raster", extended[-1, "Cell.area"], sep = "")
    output <- list(threshold = threshold,
                   extent.stand = extended[1, "Extent"],
                   occupancy.stand = extended,
                   occupancy.orig = original,
                   atlas.raster.stand = atlas_thresh,
                   scaled.rasters = scaled.rasters)
  }
  class(output) <- "upgrain"
  return(output)
}