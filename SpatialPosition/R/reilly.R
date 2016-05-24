#' @title Reilly Catchment Areas
#' @name reilly
#' @description This function computes the catchment areas as defined by W.J. Reilly (1931).
#' @param knownpts sp object (SpatialPointsDataFrame or SpatialPolygonsDataFrame); 
#' this is the set of known observations to estimate the catchment areas from.
#' @param unknownpts sp object (SpatialPointsDataFrame or SpatialPolygonsDataFrame); 
#' this is the set of unknown units for which the function computes the estimates. 
#' Not used when \code{resolution} is set up. (optional)
#' @param matdist matrix; a distance matrix. Row names match the first 
#' column of the \code{knownpts} object dataframe. Column names match the first column 
#' of the \code{unknownpts} object dataframe. (optional)
#' @param varname character; name of the variable in the \code{knownpts} dataframe from which values are computed.
#' Quantitative variable with no negative values. 
#' @param typefct character; spatial interaction function. Options are "pareto" 
#' (means power law) or "exponential".
#' If "pareto" the interaction is defined as: (1 + alpha * mDistance) ^ (-beta).
#' If "exponential" the interaction is defined as: 
#' exp(- alpha * mDistance ^ beta).
#' The alpha parameter is computed from parameters given by the user 
#' (\code{beta} and \code{span}).
#' @param span numeric; distance where the density of probability of the spatial 
#' interaction function equals 0.5.
#' @param beta numeric; impedance factor for the spatial interaction function.  
#' @param resolution numeric; resolution of the output SpatialPointsDataFrame
#'  (in map units). 
#' @param longlat logical; euclidean distance (FALSE, default) or Great Circle distance (TRUE).
#' If TRUE inputs are expected in the WGS84 reference system.
#' @param mask sp object; the spatial extent of this object is used to 
#' create the regularly spaced SpatialPointsDataFrame output. (optional)
#' @details If \code{unknownpts} is NULL then \code{resolution} must be used. 
#' @return SpatialPointsDataFrame with the computed catchment areas in a new field nammed \code{OUTPUT}.
#' Values match the row names of \code{knownpts}
#' @seealso \link{reilly}, \link{rasterReilly}, \link{plotReilly}, \link{CreateGrid}, \link{CreateDistMatrix}.
#' @examples 
#' # Create a SpatialPointsDataFrame grid of spatMask extent and 200 meters 
#' # resolution
#' data(spatData)
#' mygrid <- CreateGrid(w = spatMask, resolution = 200)
#' # Create a distance matrix between known points (spatPts) and mygrid
#' mymat <- CreateDistMatrix(knownpts = spatPts, unknownpts = mygrid, 
#'                           longlat = FALSE)
#' # Compute Reilly catchment areas from known points (spatPts) on a given 
#' # grid (mygrid) using a given distance matrix (mymat)
#' myreilly2 <- reilly(knownpts = spatPts, unknownpts = mygrid, 
#'                matdist = mymat, varname = "Capacite", 
#'                typefct = "exponential", span = 1250, 
#'                beta = 3, longlat = FALSE, mask = spatMask)
#' row.names(spatPts) <- spatPts$CodHop
#' # Compute Reilly catchment areas from known points (spatPts) on a 
#' # grid defined by its resolution
#' myreilly <- reilly(knownpts = spatPts, varname = "Capacite", 
#'                 typefct = "exponential", span = 1250, beta = 3, 
#'                 resolution = 200, longlat = FALSE, mask = spatMask)
#' # The function output a SpatialPointsDataFrame
#' class(myreilly)
#' # The OUTPUT field values match knownpts row names
#' head(unique(myreilly$OUTPUT))
#' @references REILLY, W. J. (1931) The law of retail gravitation, W. J. Reilly, New York.
#' @import sp
#' @import raster
#' @export
reilly <- function(knownpts,
                   unknownpts = NULL,
                   matdist = NULL,
                   varname,
                   typefct = "exponential", 
                   span,
                   beta,
                   resolution = 2000,
                   longlat = FALSE, 
                   mask = NULL)
{
  TestSp(knownpts)
  if (!is.null(unknownpts)){  
    TestSp(unknownpts)
    if(identicalCRS(knownpts,unknownpts) == FALSE){
      stop(paste("Inputs (",quote(knownpts), " and ",quote(unknownpts),
                 ") do not use the same projection", sep = ""),call. = FALSE)
    }
    if (!is.null(matdist)){
      matdist <- UseDistMatrix(matdist =matdist, knownpts = knownpts, 
                               unknownpts =  unknownpts) 
    }else{
      matdist <- CreateDistMatrix(knownpts = knownpts, unknownpts = unknownpts, 
                                  longlat = longlat) 
    }
  } else {
    unknownpts <- CreateGrid(w = if(is.null(mask)){knownpts} else {mask}, 
                             resolution = resolution) 
    matdist <- CreateDistMatrix(knownpts = knownpts, unknownpts = unknownpts, 
                                longlat = longlat) 
  }
  
  
  matdens <- ComputeInteractDensity(matdist = matdist, typefct = typefct,
                                    beta = beta, span = span)
  
  matopport <- ComputeOpportunity(knownpts = knownpts, matdens = matdens, 
                                  varname = varname)
  
  unknownpts <- ComputeReilly(unknownpts = unknownpts, 
                              matopport = matopport)
  
  return(unknownpts)
}

#' @title Create a Raster from a Reilly SpatialPointsDataFrame
#' @name rasterReilly
#' @description This function creates a raster from a regularly spaced 
#' Reilly SpatialPointsDataFrame (output of the \code{\link{reilly}} function). 
#' @param x sp object (SpatialPointsDataFrame); output of the \code{reilly} function.
#' @param mask sp object (SpatialPolygonsDataFrame); this object is used to clip the raster. (optional)
#' @return Raster of catchment areas values.
#' The raster uses a RAT (\code{\link{ratify}}) that contains the 
#' correspondance between raster values and catchement areas values. Use \code{
#' unique(levels(rasterName)[[1]])} to see the correpondance table.
#' @seealso \link{reilly}, \link{rasterReilly}, \link{plotReilly}, \link{CreateGrid}, \link{CreateDistMatrix}.
#' @examples 
#' data(spatData)
#' row.names(spatPts) <- spatPts$CodHop
#' # Compute Reilly catchment areas from known points (spatPts) on a
#' # grid defined by its resolution
#' myreilly <- reilly(knownpts = spatPts, varname = "Capacite",
#'                    typefct = "exponential", span = 750, beta = 2,
#'                    resolution = 50, longlat = FALSE, mask = spatMask)
#' # Create a raster of reilly values
#' myreillyraster <- rasterReilly(x = myreilly, mask = spatMask)
#' plot(myreillyraster, col = rainbow(18))
#' # Correspondance between raster values and reilly areas
#' head(unique(levels(myreillyraster)[[1]]))
#' @import sp
#' @import raster
#' @export
rasterReilly <- function(x ,mask = NULL){
  gridded(x) <- TRUE
  r <- raster(x)
  x$OUTPUT2 <- as.factor(x$OUTPUT)
  levels(x$OUTPUT2) <- 1:length(levels(x$OUTPUT2) )
  x$OUTPUT2 <- as.numeric(x$OUTPUT2)
  rasterx <- rasterize(x, r, field = 'OUTPUT2')
  if(!is.null(mask)){
    TestSp(mask)
    rasterx <- mask(rasterx, mask = mask)
  }
  ratify(rasterx)
  levels(rasterx) <- data.frame(ID = x$OUTPUT2, idarea = x$OUTPUT)
  return(rasterx)
}

#' @title Plot a Reilly Raster
#' @name plotReilly
#' @description This function plots the raster produced by the \code{\link{rasterReilly}} function.
#' @param x raster; output of the \code{\link{rasterReilly}} function.
#' @param add logical; if TRUE the raster is added to the current plot, if FALSE the raster is displayed in a new plot.
#' @param col function; color ramp function, such as \code{\link{colorRampPalette}}.
#' @details Display the raster nicely.
#' @seealso \link{reilly}, \link{rasterReilly}, \link{plotReilly}, \link{CreateGrid}, \link{CreateDistMatrix}.
#' @examples 
#' data(spatData)
#' row.names(spatPts) <- spatPts$CodHop
#' # Compute Reilly catchment areas from known points (spatPts) on a
#' # grid defined by its resolution
#' myreilly <- reilly(knownpts = spatPts, varname = "Capacite",
#'                    typefct = "exponential", span = 750, beta = 2,
#'                    resolution = 50, longlat = FALSE, mask = spatMask)
#' # Create a raster of reilly values
#' myreillyraster <- rasterReilly(x = myreilly, mask = spatMask)
#' # Plot the raster nicely
#' plotReilly(x = myreillyraster)
#' @import sp
#' @import raster
#' @importFrom grDevices rainbow
#' @export
plotReilly <- function(x, add = FALSE, 
                       col = rainbow){
  nclass <- nrow(unique(levels(x)[[1]]))
  colorReilly <- col(n = nclass)
  plot(x, legend = FALSE, axes = FALSE,
       box = FALSE, col = colorReilly,  add = add)
}