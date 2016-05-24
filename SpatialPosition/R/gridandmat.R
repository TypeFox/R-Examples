#' @title Create a Regularly Spaced SpatialPointsDataFrame
#' @name CreateGrid
#' @description This function creates a regular grid of SpatialPointsDataFrame 
#' from the extent of a given sp object and a given resolution.
#' @param w sp object; the spatial extent of this object is used to 
#' create the regular SpatialPointsDataFrame.
#' @param resolution numeric; resolution of the grid (in map units). 
#' @return The output of the function is a SpatialPointsDataFrame of regularly
#' spaced points with the same extent as \code{w}. 
#' @seealso \link{CreateDistMatrix}.
#' @examples 
#' # Create a SpatialPointsDataFrame grid of spatMask extent and 200 meters 
#' # resolution
#' data(spatData)
#' mygrid <- CreateGrid(w = spatMask, resolution = 200)
#' plot(mygrid, cex = 0.1, pch = ".")
#' plot(spatMask, border="red", lwd = 2, add = TRUE)
#' @import sp
#' @export
CreateGrid <- function (w, resolution)
{
  TestSp(w)
  boundingBox <- bbox(w)
  if(is.null(resolution)){
    resolution <- sqrt(((boundingBox[1,2] - boundingBox[1,1]) * 
                          (boundingBox[2,2] - boundingBox[2,1]))/4000)
  }
  rounder <- boundingBox %% resolution
  boundingBox[,1] <- boundingBox[,1] - rounder[,1]
  boundingBox[,2] <- boundingBox[,2] + resolution - rounder[,2]
  boxCoordX <- seq(from = boundingBox[1,1] - resolution*10, 
                   to = boundingBox[1,2]+resolution*10, 
                   by = resolution)
  boxCoordY <- seq(from = boundingBox[2,1] - resolution * 10, 
                   to = boundingBox[2,2] + resolution*10, 
                   by = resolution)
  spatGrid <- expand.grid(boxCoordX, boxCoordY)
  idSeq <- seq(1, nrow(spatGrid), 1)
  spatGrid <- data.frame(ID = idSeq, 
                         COORDX = spatGrid[, 1], 
                         COORDY = spatGrid[, 2])
  
  spatGrid <- SpatialPointsDataFrame(coords = spatGrid[ , c(2, 3)], 
                                     data = spatGrid, 
                                     proj4string = CRS(proj4string(w)))
  return(spatGrid)
}


#' @title Create a Distance Matrix Between Two Sp Objects
#' @name CreateDistMatrix
#' @description This function creates a distance matrix between two 
#' sp objects (SpatialPointsDataFrame or SpatialPolygonsDataFrame).
#' @param knownpts sp object; rows of the distance matrix.
#' @param unknownpts sp object; columns of the distance matrix.
#' @param longlat logical; euclidean distance (FALSE, default) or Great Circle distance (TRUE).
#' @param bypassctrl logical; bypass the distance matrix size control (see Details).
#' @details The function returns a full matrix of distances in the metric of the 
#' points if \code{longlat} is FALSE, or in kilometers if \code{longlat} is TRUE. This is a wrapper
#' for the \code{\link{spDists}} function. 
#' 
#' If the matrix to compute is too large (more than 100,000,000 cells or more than 10,000,000 origins or destinations) 
#' the function sends a confirmation message to warn users about the amount of RAM mobilized. 
#' Use \code{bypassctrl} = TRUE to skip this control.
#' @return A distance matrix, row names are \code{knownpts} row names, column names are \code{unknownpts} row names.
#' @seealso \link{CreateGrid}.
#' @examples 
#' # Create a SpatialPointsDataFrame grid of spatMask extent and 200 meters 
#' # resolution
#' data(spatData)
#' mygrid <- CreateGrid(w = spatMask, resolution = 200)
#' # Create a distance matrix between known spatPts and mygrid
#' mymat <- CreateDistMatrix(knownpts = spatPts, unknownpts = mygrid, 
#'                           longlat = FALSE, bypassctrl = FALSE)
#' mymat[1:5,1:5]
#' nrow(spatPts)
#' nrow(mygrid)
#' dim(mymat)
#' @import sp
#' @export
CreateDistMatrix  <- function(knownpts, 
                              unknownpts, 
                              longlat = FALSE, 
                              bypassctrl = FALSE)
{
  TestSp(knownpts)
  TestSp(unknownpts)
  if(identicalCRS(knownpts,unknownpts) == FALSE){
    stop(paste("Inputs (",quote(knownpts), " and ",quote(unknownpts),
               ") do not use the same projection", sep = ""),call. = FALSE)
  }
  if (bypassctrl == FALSE){
    nk <- nrow(knownpts)
    nu <- nrow(unknownpts)
    if(nk * nu > 100000000 | nu > 10000000 | nk > 10000000){
      if (interactive()){
        cat("Do you really want to this distance matrix (from", nk , 
            "known points to", nu,"estimated values) ? \n 
            (It seems to be a heavy computation.) [y/n]" )
        z <- readLines(con = stdin(), n = 1) 
        while (!z %in% c("n","y")){
          cat ("Enter y or n")
          z <- readLines(con = stdin(), n = 1)  
        }
        if (z == "y"){
          cat("Ok, YOLO!")
        } else {
          stop("Computation aborted. Matrix would probably be too big.",
               call. = F)
        }
      } else {
        stop("Computation aborted. Matrix would probably be too big.", 
             call. = F)
      }
    }
  }
  matDist <- spDists(x = knownpts, y = unknownpts, longlat = longlat)
  dimnames(matDist) <- list(row.names(knownpts), row.names(unknownpts))
  return(round(matDist, digits = 8))
}
