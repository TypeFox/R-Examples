#' @title Stewart Potentials
#' @name stewart
#' @description This function computes the potentials as defined by J.Q. Stewart (1942).
#' @param knownpts sp object (SpatialPointsDataFrame or SpatialPolygonsDataFrame);
#' this is the set of known observations to estimate the potentials from.
#' @param unknownpts sp object (SpatialPointsDataFrame or SpatialPolygonsDataFrame); 
#' this is the set of unknown units for which the function computes the estimates. 
#' Not used when \code{resolution} is set up. (optional)
#' @param matdist matrix; a distance matrix. Row names match the first 
#' column of the \code{knownpts} object dataframe. Column names match the first column 
#' of the \code{unknownpts} object dataframe. (optional)
#' @param varname character; name of the variable in the \code{knownpts} dataframe from which potentials are computed.
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
#' @return SpatialPointsDataFrame with the computed potentials in a new field nammed \code{OUTPUT}
#' @seealso \link{rasterStewart}, \link{plotStewart}, \link{quickStewart},
#' \link{rasterToContourPoly}, \link{CreateGrid}, \link{CreateDistMatrix}.
#' @examples 
#' # Create a SpatialPointsDataFrame grid of spatMask extent and 200 meters 
#' # resolution
#' data(spatData)
#' mygrid <- CreateGrid(w = spatMask, resolution = 200)
#' # Create a distance matrix between known points (spatPts) and mygrid
#' mymat <- CreateDistMatrix(knownpts = spatPts, unknownpts = mygrid, 
#'                           longlat = FALSE)
#' # Compute Stewart potentials from known points (spatPts) on a given 
#' # grid (mygrid) using a given distance matrix (mymat)
#' mystewart <- stewart(knownpts = spatPts, unknownpts = mygrid, 
#'                      matdist = mymat, varname = "Capacite", 
#'                      typefct = "exponential", span = 1250, 
#'                      beta = 3, longlat = FALSE, mask = spatMask)
#' # Compute Stewart potentials from known points (spatPts) on a 
#' # grid defined by its resolution
#' mystewart2 <- stewart(knownpts = spatPts, varname = "Capacite", 
#'                       typefct = "exponential", span = 1250, beta = 3, 
#'                       resolution = 200, longlat = FALSE, mask = spatMask)
#' # The two methods have the same result
#' identical(mystewart, mystewart2)
#' # the function output a SpatialPointsDataFrame
#' class(mystewart)
#' # Computed values
#' summary(mystewart$OUTPUT)
#' @references 
#' STEWART J.Q. (1942) "Measure of the influence of a population at a distance", Sociometry, 5(1): 63-71.  
#' @import sp
#' @import raster
#' @export
stewart <- function(knownpts,
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
  
  unknownpts <- ComputePotentials(unknownpts = unknownpts, 
                                  matopport = matopport)
  
  return(unknownpts)
}

#' @title Create a Raster from a Stewart SpatialPointsDataFrame
#' @name rasterStewart
#' @description This function creates a raster from a regularly spaced 
#' Stewart SpatialPointsDataFrame (output of the \code{\link{stewart}} function). 
#' @param x sp object (SpatialPointsDataFrame); output of the \code{stewart} function.
#' @param mask sp object (SpatialPolygonsDataFrame); this object is used to clip the raster. (optional)
#' @return Raster of potential values.
#' @seealso \link{stewart}, \link{quickStewart}, \link{plotStewart}, \link{rasterToContourPoly}, \link{CreateGrid}, \link{CreateDistMatrix}.
#' @examples
#' data(spatData)
#' # Compute Stewart potentials from known points (spatPts) on a
#' # grid defined by its resolution
#' mystewart <- stewart(knownpts = spatPts, varname = "Capacite",
#'                      typefct = "exponential", span = 1000, beta = 3,
#'                      resolution = 50, longlat = FALSE, mask = spatMask)
#' # Create a raster of potentials values
#' mystewartraster <- rasterStewart(x = mystewart, mask = spatMask)
#' plot(mystewartraster)
#' @import sp
#' @import raster
#' @export
rasterStewart <- function(x, mask = NULL){
  gridded(x) <- TRUE
  r <- raster(x)
  rasterx <- rasterize(x, r, field = 'OUTPUT')
  if(!is.null(mask)){
    TestSp(mask)
    rasterx <- mask(rasterx, mask = mask)
  }
  return(rasterx)
}




#' @title Plot a Stewart Raster
#' @name plotStewart
#' @description This function plots the raster produced by the \code{\link{rasterStewart}} function.
#' @param x raster; output of the \code{\link{rasterStewart}} function.
#' @param add logical; if TRUE the raster is added to the current plot, if FALSE the raster is displayed in a new plot.
#' @param breaks numeric; vector of break values to map. If used, 
#' this parameter overrides \code{typec} and \code{nclass} parameters 
#' @param typec character; either "equal" or "quantile", how to discretize the values.
#' @param nclass numeric (integer), number of classes.
#' @param legend.rnd numeric (integer); number of digits used to round the values displayed in the legend.
#' @param col function; color ramp function, such as \code{\link{colorRampPalette}}.
#' @return Display the raster nicely and return the list of break values (invisible).
#' @seealso \link{stewart}, \link{rasterStewart}, \link{quickStewart}, \link{rasterToContourPoly}, \link{CreateGrid}, \link{CreateDistMatrix}.
#' @examples 
#' data(spatData)
#' # Compute Stewart potentials from known points (spatPts) on a
#' # grid defined by its resolution
#' mystewart <- stewart(knownpts = spatPts, varname = "Capacite",
#'                      typefct = "exponential", span = 1000, beta = 3,
#'                      resolution = 50, longlat = FALSE, mask = spatMask)
#' # Create a raster of potentials values
#' mystewartraster <- rasterStewart(x = mystewart, mask = spatMask)
#' # Plot stewart potentials nicely
#' plotStewart(x = mystewartraster, add = FALSE, nclass = 5)
#' # Can be used to obtain break values
#' break.values <- plotStewart(x = mystewartraster, add = FALSE, nclass = 5)
#' break.values
#' @import sp
#' @import raster
#' @importFrom grDevices colorRampPalette
#' @export
plotStewart <- function(x, add = FALSE, 
                        breaks = NULL, typec = "equal", 
                        nclass = 5, legend.rnd = 0, 
                        col =  colorRampPalette(c("#FEA3A3","#980000"))){
  if (!is.null(breaks)){
    bks <- unique(breaks[order(breaks)])
  } else if (typec == "equal"){
    bks <- seq(from = cellStats(x, min), 
               to = cellStats(x, max), length.out = nclass+1)
  } else if (typec == "quantile"){
    bks <- quantile (x, probs = seq(0,1, by = 1/nclass))
  } else {
    stop('Enter a proper discretisation type: "equal" or "quantile"')
  }
  bks <- unique(bks)
  col <- col(length(bks)-1)
  plot(x, breaks = bks, legend = FALSE, axes = FALSE,
       box = FALSE, col = col,  add = add)
  
  nbks <- round(bks,legend.rnd)
  leglab <- rep(NA, (length(nbks)-1))
  for(i in 1:(length(nbks)-1)){
    leglab[i] <- paste("[", nbks[i], " - ", nbks[i+1],"[" ,sep="")
  }
  leglab[i] <- paste( substr(leglab[i],1, nchar(leglab[i])-1), "]", sep="")
  
  graphics::legend(x='topright', legend = rev(leglab), 
                   xpd=T,inset=c(-0.2,0), 
                   fill = rev(col), cex = 0.7, 
                   plot = TRUE, bty = "n", 
                   title = "Potentials")
  
  
  
  #   plot(x, legend.only=TRUE, col = col, 
  #        breaks=round(bks,legend.rnd))
  
  return(invisible(bks))
}

