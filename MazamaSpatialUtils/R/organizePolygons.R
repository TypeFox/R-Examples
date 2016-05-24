#' @keywords internal
#' @export
#' @title Organize Ungrouped Polygons
#' @param dataset object of class SpatialPolygonsDataFrame
#' @param uniqueID unique identifier to determine which values are duplicated
#' @param sumColumns vector of column names to be summed
#' @description Determines if the SpatialPolygonsDataFrame is grouped. If ungrouped, 
#' function will group duplicated values based on the provided unique identifier. 
#' 
#' If sumColumns is NULL and there are multiple rows that aren't duplicated but have
#' the same 'uniqueID', the original SpatialPolygonsDataFrame will be returned.
#' @return SpatialPolygonsDataFrame composed of grouped polygons.
#' @examples SPDF <- organizePolygons(SimpleTimezones, 'timezone', NULL)
#' @examples SPDF <- organizePolygons(USGS_HUC_8, 'huc', 'area')
organizePolygons <- function(dataset, uniqueID, sumColumns=NULL) {
  
  # Test if the unique identifier is a character string
  if (!is.character(uniqueID)) {
    stop(paste0('The uniqueID, "',uniqueID,'" must be a character string.'), call.=FALSE)
  }
  
  # Test if the dataset is of class SpatialPolygonsDataFrame
  if (!class(dataset)=='SpatialPolygonsDataFrame') {
    stop(paste0(dataset,' not found. Please use loadSpatialData().'), call.=FALSE)    
  }
    
  # Test if the dataframe contains grouped polygons, return the dataframe if true
  if ( !any(duplicated(dataset@data[,uniqueID])) ) {
    return(dataset)
  }
  
  # Determine which values are duplicated and create an updated dataframe
  if (is.null(sumColumns)) {
    dupMask <- duplicated(dataset@data)    
  } else {
    dupMask <- duplicated(dataset@data[,uniqueID])
  }
  nonDups <- dataset[!dupMask,]
  
  # Test if there are any rows that have non-duplicated data but no columns specified for summation
  if (any(duplicated(nonDups@data[,uniqueID])) && is.null(sumColumns)) {
    message(paste0('There are duplicated ',uniqueID,' rows with different values. ',
                   'Please specify columns to be summed. Returning original dataframe.'))
    return(dataset)
  }
  
  # Create an empty list to store each list Polygons
  srl <- list()
  
  # Group polygons based off the unique identifier
  for (i in seq_along(nonDups)) {
    ###print(i)
    x <- nonDups@data[,uniqueID][i]
    allX <- which(dataset@data[,uniqueID] == x)
    
    # Create an emply list to store the Polygons corresponding 'x'
    newPolygons <- list()
    for(index in seq_along(allX)) {
      newPolygons[[index]] <- dataset@polygons[[ allX[index] ]]@Polygons[[1]]
    }
    
    # Create an object of class Polygons
    polys <- sp::Polygons(newPolygons, x)
    srl[[i]] <- polys
    
    # If a vector of column names is given, sum up those columns and replace the old row
    if (!is.null(sumColumns)) {
      for (j in seq_along(sumColumns)) {
        nonDups@data[ i, sumColumns[j] ] <- sum(dataset@data[ allX, sumColumns[j] ])
      }
    }
  }
  
  # Create an object SpatialPolygons
  proj4 <- dataset@proj4string
  SP <- sp::SpatialPolygons(srl, proj4string=proj4)
  
  # Build a SpatialPolygonsDataFrame from the dataframe and SpatialPolygons
  rownames(nonDups@data) <- nonDups@data[,uniqueID]
  SPDF <- sp::SpatialPolygonsDataFrame(SP, nonDups@data)
  
  return(SPDF)
}

