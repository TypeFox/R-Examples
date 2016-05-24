#' @title Raster Entropy
#' @description Calculates entropy on integer raster (i.e., 8 bit 0-255)  
#'                                                                       
#' @param x object of class raster (requires integer raster)  
#' @param d Size of matrix (window)
#' @param filename Raster file written to disk
#' @param ... Optional arguments passed to writeRaster or dataType              
#'  
#' @return raster class object or specified format raster written to disk                
#'
#' @note
#' Entropy calculated as: H = -sum(Pi*ln(Pi))
#'   where; Pi, Proportion of one value to total values Pi=n(p)/m and m, Number of unique values
#'   Expected range: 0 to log(m)
#'     H=0 if window contains the same value in all cells.
#'     H increases with the number of different values in the window.
#'
#' Maximum entropy is reached when all values are different, same as log(m)
#'   max.ent <- function(x) { log( length( unique(x) ) ) }
#'
#' @note Depends: raster
#'
#' @references
#' Fuchs M., Hoffmann R., Schwonke F. (2008) Change Detection with GRASS GIS - Comparison of images taken by different sensor. On line at: http://geoinformatics.fsv.cvut.cz/gwiki/Change_Detection_with_GRASS_GIS_-_Comparison_of_images_taken_by_different_sensors
#'
#' @examples 
#' require(raster)
#'   r <- raster(ncols=100, nrows=100)
#'     r[] <- round(runif(ncell(r), 1,8), digits=0)
#'
#' rEnt <- raster.entropy(r, d=5)
#'   opar <- par  
#'     par(mfcol=c(2,1))
#'       plot(r)
#'         plot(rEnt)
#'   par(opar)
#'
#' @export  
raster.entropy <- function(x, d = 5, filename = FALSE, ...) {
    if (!inherits(x, "RasterLayer")) 
        stop("MUST BE RasterLayer OBJECT")
    entropy <- function(x) {
        p <- vector()
        if (length(unique(x)) <= 1) {
            return(0)
        }
        nv <- length(unique(x))
        for (i in unique(x)) {
            p <- append(p, (length(x[x == i])/nv * log(length(x[x == i])/nv)))
        }
        return(-sum(p))
    }
    if (filename != FALSE) {
        raster::focal(x, w = d, fun = entropy, filename = filename, ...)
        print(paste("RASTER WRITTEN TO", filename, sep = ": "))
    } else {
        return(raster::focal(x, w = matrix(1, nrow = d, ncol = d), fun = entropy))
    }
} 
