#!/usr/bin/env Rscript

#' rescaleLayer
#'
#' rescale the RasterLayer values to min:0 max:1
#'
#' @param raster.layer an object of RasterLayer class
#' @return an object of RasterLayer that rescaled.
rescaleLayer <- function(raster.layer) {
    if (!(class(raster.layer) %in% "RasterLayer")) {
        stop("raster.layer is not a RasterLayer objectect!")
    }
    min.value <- cellStats(raster.layer, min)
    if (min.value < 0) {
        raster.layer <- raster.layer + (0-min.value)
        max.value <- cellStats(raster.layer, max)
        raster.layer <- raster.layer/max.value
        return(raster.layer)
    } else {
        max.value <- cellStats(raster.layer, max)
        raster.layer <- raster.layer/max.value
        return(raster.layer)
    }
}

#' rescaleStack
#'
#' rescale the RasterStack values to min:0 max:1
#'
#' @param raster.stack an object of RasterStack class
#' @return an object of RasterStack that rescaled.
rescaleStack <- function(raster.stack) {
    if (!(class(raster.stack) %in% "RasterStack")) {
        stop("env.stack is not a RasterStack object!")
    }
    raster.name <- names(raster.stack)
    raster.stack.list <- lapply(X=raster.name, FUN=function(name, raster.stack) {return(rescaleLayer(raster.stack[[name]]))}, raster.stack)
    result.stack <- stack(raster.stack.list)
    names(result.stack) <- raster.name
    return(result.stack)
}

#' rescale
#'
#' rescale the RasterStack or RasterLayer values to min:0 max:1
#'
#' @param raster.object an object of RasterStack or RasterLayer class
#' @return an object of RasterStack or RasterLayer that rescaled.
#' @export
rescale <- function(raster.object) {
    if (!(class(raster.object) %in% c("RasterLayer", "RasterStack"))) {
        stop("raster.object should be a RasterLayer or RasterStack object!")
    }
    if (class(raster.object) %in% "RasterLayer") {
        raster.object <- rescaleLayer(raster.object)
    } else {
        raster.object <- rescaleStack(raster.object)
    }
    return(raster.object)
}
