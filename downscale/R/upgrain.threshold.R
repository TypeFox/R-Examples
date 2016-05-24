################################################################################
#
# upgrain.threshold.R 
# Version 1.1
# 05/10/2015
#
# Updates:
#   05/10/2015: Bug on number of scales fixed
#               'All_Presences' changed to 'All_occurrences'
#
# Explores the NoData threshold selection of cell selection for upgraining. A 
# threshold of 0 means that all cells at the largest grain with some Sampled 
# are included (All Sampled). A threshold of 1 means only cells at the largest 
# grain that are entirely within the Sampled are included (Sampled only). Two 
# other threshold options are also suggested based upon the trade-off of 
# including NoData cells as absence and losing presence records: 
# All Occurences: the threshold that preserves all presence records
# Gain Equals Loss: the threshold where the upgrained extent = the original 
# extent
#
# args:
#   atlas.data = either a raster file of presence-absence atlas data, or a data 
#                frame of sampled cells with longitude, latitude and 
#                presence-absence.
#   cell.width = if data is a data frame, the cell widths of sampled cells.
#   scales = the number of cells to upgrain. Upgraining will happen by factors 
#            of 2 - ie if scales = 3, the atlas data will be aggregated in 2x2 
#            cells, 4x4 cells and 8x8 cells.
#   thresholds = a vector of thresholds between and including 0 and 1 for the
#                quantity of Unsampled NoData cells that can be included.
#
# outputs:
#   Two plotting windows. One with four plots relavent to exploring the 
#   trade-offs. The second plots the maps of the four default threshold 
#   selections.
#   A list with two objects;
#     Thresholds: the threshold values for the four default threshold selections
#     Data: Data frame with all the relevant data for exploring the trade-offs.
#
################################################################################

upgrain.threshold <- function(atlas.data,
                              cell.width  = NULL,
                              scales,
                              thresholds = seq(0, 1, 0.01)) {
  ### Error checking: thresholds are between and include 0 and 1
  if(min(thresholds) != 0){
    stop("Threshold values must be between and include 0 and 1")
  }
  if(max(thresholds) != 1){
    stop("Threshold values must be between and include 0 and 1")
  }
  
  ### Error checking: if data frame needs cell width
  if(is.data.frame(atlas.data) == "TRUE") {
    if(is.null(cell.width)) {
      stop("If data is data.frame cell.width is required")
    }
  }
  
  ### data manipulation
  if(is.data.frame(atlas.data) == "TRUE") {
    ### error checking: format of data frame
    if(ncol(atlas.data) != 3) {
      stop("Input data frame must contain three columns named lon, lat, 
           presence in that order")
    }
    ### error checking: column names
    if(sum(names(atlas.data) != c("lon", "lat", "presence")) > 0) {
      stop("Input data frame must contain three columns named lon, lat, 
           presence in that order")
    }
    
    longitude <- atlas.data[, "lon"]
    latitude <- atlas.data[, "lat"]
    presence <- atlas.data[, "presence"]
    presence[presence > 0] <- 1
    shapefile <- sp::SpatialPointsDataFrame(coords = data.frame(lon = longitude,
                                                                lat = latitude),
                                            data = data.frame(presence = 
                                                                presence)) 
    atlas_raster <- raster::raster(ymn = min(latitude) - (cell.width / 2),
                                   ymx = max(latitude) + (cell.width / 2),
                                   xmn = min(longitude) - (cell.width / 2), 
                                   xmx = max(longitude) + (cell.width / 2),
                                   resolution = cell.width)
    atlas_raster <- raster::rasterize(shapefile, atlas_raster)
    atlas_raster <- raster::dropLayer(atlas_raster, 1)
  }
  
  if(class(atlas.data) == "RasterLayer") {
    atlas_raster <- atlas.data
    atlas_raster@data@values[atlas_raster@data@values > 0] <- 1
  }
  
  ### largest scale (all other layers are extended to equal this raster)
  max_raster <- raster::aggregate(atlas_raster, (2 ^ scales), fun = max)
  
  ########################################################
  ### create boundary raster
  boundary_raster <- atlas_raster
  boundary_raster@data@values[!is.na(boundary_raster@data@values)] <- 1
  boundary_raster <- raster::aggregate(boundary_raster, (2 ^ scales), fun = sum)  
  boundary_raster@data@values <- boundary_raster@data@values / ((2 ^scales )^ 2)
  boundary_poly <- raster::rasterToPolygons(boundary_raster, dissolve = TRUE)
  
  ##################################################
  ### Loop through thresholds
  thresholds <- thresholds
  land <- data.frame(Threshold = thresholds, 
                     SampledExcluded = rep(NA, length(thresholds)),
                     SampledIncluded = NA,
                     UnsampledAdded = NA, 
                     Extent = NA,
                     OccurrencesExcluded = NA)
  
  atlas_raster_extend <- raster::extend(atlas_raster, 
                                        raster::extent(max_raster))
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
  ### Possible Thresholds
  Gain_loss.thresh <- thresholds[which.min(abs(land[, "Extent"] - 
                                                 sum(atlas_boundary@data@values,
                                                     na.rm = TRUE)))]
  Presence.thresh <- thresholds[max(which(land[, "OccurrencesExcluded"] == 0))]
  
  ### Plotting
  par.original <- par()
  par.original <- list(mfrow = par.original$mfrow, mar = par.original$mar)
  par(mfrow=c(2,2), mar = c(5,5,1,1))
  
  # Plot 1: the total extent
  plot(land[, "Threshold"],
       land[, "Extent"],
       type = "l",
       ylab="Extent (total number of cells)",
       xlab="Threshold",
       ylim = c(0, max(land[, "Extent"], na.rm = TRUE)))
  abline(h = sum(atlas_boundary@data@values == 1, na.rm = TRUE),
         col = "Grey")
  lines(c(Gain_loss.thresh, Gain_loss.thresh),
        c(0, land[thresholds == Gain_loss.thresh, "Extent"]),
        col = "red", lty = 2)
  lines(c(Presence.thresh, Presence.thresh),
        c(0, land[thresholds == Presence.thresh, "Extent"]),
        col = "blue", lty = 2)
  legend("topright",
         legend = c("Original extent",
                    "Gain equals loss threshold",
                    "All occurrences threshold"),
         lty = c(1,2,2), bty = "n", col = c("Grey", "Red", "Blue"))
  
  # Plot 2: the number of Sampled cells incorrectly identified
  plot(land[, "Threshold"],
       land[, "SampledExcluded"],
       type = "l",
       ylab="Number of cells",
       xlab="Threshold",
       ylim = c(0, max(c(land[, "SampledExcluded"],
                         land[, "UnsampledAdded"]), na.rm = TRUE)))
  lines(land[, "Threshold"], land[, "UnsampledAdded"], lty = 2)
  legend("top",
         legend = c("Number of sampled cells excluded",
                    "Number of unsampled cells added"),
         lty = 1:2, bty = "n")
  lines(c(Gain_loss.thresh, Gain_loss.thresh),
        c(0, land[thresholds == Gain_loss.thresh, "SampledExcluded"]),
        col = "red", lty = 2)
  lines(c(Presence.thresh, Presence.thresh),
        c(0, land[thresholds == Presence.thresh, "UnsampledAdded"]),
        col = "blue", lty = 2)
  
  # Plot 3: the number of original land cells retained
  plot(land[, "Threshold"],
       land[, "SampledIncluded"] / max(land[, "SampledIncluded"],
                                        na.rm = TRUE),
       type = "l",
       ylab="Prop. of sampled cells retained",
       xlab="Threshold",
       ylim = c(0, 1))
  lines(c(Gain_loss.thresh, Gain_loss.thresh),
        c(0, (land[thresholds == Gain_loss.thresh, "SampledIncluded"] /
                max(land[, "SampledIncluded"]))),
        col = "red", lty = 2)
  lines(c(Presence.thresh, Presence.thresh),
        c(0, (land[thresholds == Presence.thresh, "SampledIncluded"] /
                max(land[, "SampledIncluded"]))),
        col = "blue", lty = 2)
  
  # Plot 4: % of presence cells exluded
  plot(land[, "Threshold"],
       land[, "OccurrencesExcluded"],
       type = "l",
       ylab="Prop. of Occurrences excluded",
       xlab="Threshold",
       ylim = c(0, 1))
  lines(c(Gain_loss.thresh, Gain_loss.thresh),
        c(1, land[thresholds == Gain_loss.thresh, "OccurrencesExcluded"]),
        lty = 2, col = "red")
  lines(c(Presence.thresh, Presence.thresh),
        c(1, land[thresholds == Presence.thresh, "OccurrencesExcluded"]),
        lty = 2, col = "blue")
  
  ### Maps of thresholds options
  par(ask = TRUE) 
  par(mfrow=c(2,2), mar = c(5.5,1,3.5,1))
  
  thresh.selection <- c(0, Presence.thresh, Gain_loss.thresh, 1)
  selection <- c("All Sampled",
                 "All Occurrences",
                 "Gain Equals Loss",
                 "Sampled Only")
  selection <- selection[order(thresh.selection)]
  thresholds <- thresholds[order(thresh.selection)] 
  
  for(j in 1:length(thresh.selection)){
    max_raster_thresh <- max_raster
    max_raster_thresh@data@values[boundary_raster@data@values <
                                    thresh.selection[j]] <- NA
    
    max_raster_thresh_extent <- max_raster_thresh
    max_raster_thresh_extent@data@values[max_raster_thresh_extent@data@values ==
                                           0] <- 1
    boundary_poly <- raster::rasterToPolygons(max_raster_thresh_extent,
                                              dissolve = TRUE)
    
    max_raster_thresh <- raster::disaggregate(max_raster_thresh, 8, fun = max)
    
    atlas_boundary <- atlas_raster
    atlas_boundary@data@values[atlas_boundary@data@values == 0] <- 1
    atlas_boundary@data@values[is.na(atlas_boundary@data@values)] <- 0
    atlas_boundary <- raster::extend(atlas_boundary, 
                                     raster::extent(max_raster_thresh))
    
    atlas_boundary_thresh <- max_raster_thresh
    atlas_boundary_thresh@data@values[atlas_boundary_thresh@data@values >=
                                        0] <- 2
    atlas_boundary_thresh@data@values[is.na(atlas_boundary_thresh@data@values)] <- 0
    
    plot(atlas_raster_extend, axes = FALSE,
         colNA = rgb(0.5, 0.5, 0.5),
         col = c(rgb(1, 1, 1), rgb(1, 0, 0)),
         legend = FALSE)
    title(main = paste(selection[j], "\n", "Threshold = ", thresh.selection[j]),
          line = 1)
    title(sub = paste("Prop. Sampled cells retained = ",
                      round(land[land[, 1] == thresh.selection[j],
                                 "SampledIncluded"] / 
                              max(land[, "SampledIncluded"], na.rm = TRUE), 3),
                      "\n",
                      "Prop. Occurrences excluded = ", 
                      round(land[land[ ,1] == thresh.selection[j],
                                 "OccurrencesExcluded"], 2),
                      "\n",
                      "Sampled cells excuded = ", 
                      land[land[, 1] == thresh.selection[j], "SampledExcluded"],
                      "\n", 
                      "Unsampled cells added = ", 
                      land[land[, 1] == thresh.selection[j], "UnsampledAdded"],
                      sep = ""),
          line = 3.5)
    plot(atlas_boundary_thresh, 
         colNA = rgb(0.5, 0.5, 0.5),
         col = c(rgb(0.5, 0.5, 0.5), rgb(1, 0, 0)),
         alpha = 0.35,
         add = TRUE,
         legend = FALSE)
    plot(boundary_poly, add = TRUE)
  }
  par(mfrow = par.original$mfrow, mar = par.original$mar, ask = FALSE)
  
  output <- list(Thresholds = data.frame(All_Sampled = 0,
                                         All_Occurrences = Presence.thresh,
                                         Gain_Equals_Loss = Gain_loss.thresh,
                                         Sampled_Only = 1),
                 Data =  land)
  return(output)
}