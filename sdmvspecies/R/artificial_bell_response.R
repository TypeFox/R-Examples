#!/usr/bin/env Rscript

.artificialGaussianTranslate <- function(factor, normal.mean, normal.sd, rescale) {
    result <- dnorm(factor, normal.mean, normal.sd)
    if (rescale) {
        result <- sqrt(2 * pi) * normal.sd * result
    }
    return(result)
}

#' artificialBellResponse
#'
#' artificial bell response method
#'
#' This method mainly implement artificial bell response method, more detail see references.
#'
#' @param env.stack a \code{rasterStack} object that contain the environment variable
#' @param config config is a \code{list} or \code{matrix} or \code{data.frame} that contain config info, details see details part
#' @param stack stack is an option that if you want not compose them togethor (result return as a \code{rasterStack}). Default is FALSE
#' @param compose the method compose the suitability together. Default is product
#' @param rescale if \code{TRUE} each environment will rescale before compose together
#' @return \code{rasterLayer} or \code{rasterStack} if stack is set to TRUE
#' @references Varela, S., Anderson, R. P., García-Valdés, R., & Fernández-González, F. (2014). Environmental filters reduce the effects of sampling bias and improve predictions of ecological niche models. Ecography.
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
#' file.name <- files <- c("bio1.bil", "bio12.bil", "bio7.bil", "bio5.bil")
#' files <- paste(env.dir, file.name, sep="")
#' # make raster stack
#' env.stack <- stack(files)
#' # config
#' config <- list(c("bio1",150, 50), c("bio12", 2000, 500), c("bio7", 400, 100), c("bio5", 300, 100))
#' # run pick mean
#' species.raster <- artificialBellResponse(env.stack, config)
#' # plot map
#' plot(species.raster)
#' # species distribution map
#' species.distribution.raster <- species.raster > 0.2
#' # plot map
#' plot(species.distribution.raster)
artificialBellResponse <- function(env.stack, config, stack=FALSE, compose="product", rescale=TRUE) {
    # check env.stack first
    if (! (class(env.stack) %in% "RasterStack")) {
        stop("env.stack is not a RasterStack object!")
    }
    # TODO:here used mclapply but not given core.number
    species.list <- lapply(X=config, FUN=.artificialBellResponseMain, env.stack, rescale)

    species.matrix <- matrix(unlist(species.list), ncol=length(config), byrow=FALSE)
    if (! stack) {
        if (compose == "product") {
            species <- apply(species.matrix, 1, prod)
        } else if (compose == "sum") {
            species <- apply(species.matrix, 1, sum)
        }
        species.layer <- env.stack[[config[[1]][1]]]
        species.raster <- setValues(species.layer, as.vector(species))
        return(species.raster)
    } else {
        species.layer <- env.stack[[config[[1]][1]]]
        col.number <- ncol(species.matrix)
        species.stack <- stack()
        for (col.index in 1:col.number) {
            species.raster <- setValues(species.layer, as.vector(species.matrix[,col.index]))
            species.stack <- stack(species.stack, species.raster)
        }
        config.matrix <- matrix(unlist(config), byrow=TRUE, ncol=3)
        layer.names <- config.matrix[,1]
        names(species.stack) <- layer.names
        return(species.stack)
    }
}

.artificialBellResponseMain <- function(var, env.stack, rescale) {
    predictor.name <- var[1]
    normal.mean <- var[2]
    normal.mean <- as.integer(normal.mean)
    normal.sd <- var[3]
    normal.sd <- as.integer(normal.sd)

    env.layer <- env.stack[[predictor.name]]

    factor <- getValues(env.layer)

    result <- .artificialGaussianTranslate(factor, normal.mean, normal.sd, rescale)
    return(result)
}
