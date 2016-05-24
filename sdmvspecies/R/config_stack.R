#!/usr/bin/env Rscript

#' configStack
#'
#' output config layers as rasterStack
#'
#' This method will extract rasterLayer acorrding to config, then output rasterStack as result
#'
#' @param env.stack a \code{rasterStack} object that contain the environment variable
#' @param config config is a \code{list} or \code{matrix} or \code{data.frame} that contain config info, details see details part
#' @return \code{rasterStack} object
#' @encoding utf-8
#' @export
#' @examples
#' # load the sdmvspecies library
#' library("sdmvspecies")
#' library("raster")
#' # find package's location
#' package.dir <- system.file(package="sdmvspecies")
#' # let see where is our sdmvspecies is installed in
#' package.dir
#' # find env dir under the package's location
#' env.dir <- paste(package.dir, "/external/env/", sep="")
#' # let see env dir
#' env.dir
#' # get the environment raster file
#' env.files <- list.files(env.dir, pattern="*.bil$", full.names=TRUE)
#' # see the file list
#' env.files
#' # put the environment file in a raster stack,
#' # which require all the environment should have same resolution and extend
#' env.stack <- stack(env.files)
#' # let see the env.stack var
#' env.stack
#' # here let's configure the environment response function and weight
#' config <- list(c("bio1", "1", 0, 100), c("bio11", "2", -100, NULL))
#' env.stack <- configStack(env.stack, config)
#' plot(env.stack)
configStack <- function(env.stack, config) {
    # check env.stack first
    if (!(class(env.stack) %in% "RasterStack")) {
        stop("env.stack is not a RasterStack object!")
    }
    # check params and get raster name
    raster.name <- c()
    for (config.item in config) {
        if (!(config.item[1] %in% names(env.stack))) {
            stop("params must match with env.stack in layer names")
        }
        raster.name <- c(raster.name, config.item[1])
    }
    # TODO:here used mclapply but not given core.number
    raster.list <- mclapply(X=config, FUN=.extractRaster, env.stack)
    #print(species.list)

    raster.stack <- stack(raster.list)
    names(raster.stack) <- raster.name
    return(raster.stack)
}


.extractRaster <- function(var, env.stack) {
    raster.name <- var[1]
    raster.layer <- env.stack[[raster.name]]
    return(raster.layer)
}
