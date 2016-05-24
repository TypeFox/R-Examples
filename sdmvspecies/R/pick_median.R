#!/usr/bin/env Rscript


#' pickMedian
#'
#' pick median method
#'
#' This method mainly implement pick median method
#'
#' @param env.stack a \code{rasterStack} object that contain the environment variable
#' @param subset subset is a string \code{vector} that contain environment variables names which into calculate, if NULL that all var in env.stack will calculate.
#' @param stack stack is an option that if you want not compose them togethor (result return as a \code{rasterStack}). Default is FALSE
#' @return \code{rasterLayer} or \code{rasterStack} if stack is set to TRUE
#' @references Lobo, J. M., & Tognelli, M. F. (2011). Exploring the effects of quantity and location of pseudo-absences and sampling biases on the performance of distribution models with limited point occurrence data. Journal for Nature Conservation, 19(1), 1-7.
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
#' files <- list.files(path=env.dir, pattern="*.bil$", full.names=TRUE)
#' # make raster stack
#' env.stack <- stack(files)
#' # run pick mean
#' species.raster <- pickMedian(env.stack)
#' # plot map
#' plot(species.raster)
pickMedian <- function(env.stack, subset=NULL, stack=FALSE) {
    if (!(class(env.stack) %in% "RasterStack")) {
        stop("env.stack is not a RasterStack object!")
    }
    env.names <- unlist(names(env.stack))

    if (length(subset)) {
        check.result <- subset %in% env.names
        if (!all(check.result)) {
            mis.name <- subset[!check.result]
            message <- paste("Follow variable:", paste(mis.name, collapse=","), "are missing!\nAlready auto removed!", sep="")
            warning(message)
        }
        env.names <- subset[check.result]
    }

    # TODO:here used mclapply but not given core.number
    species.list <- mclapply(X=env.names, FUN=.pickMedian, env.stack)

    if (!stack) {
        species.stack <- stack(species.list)
        species.layer <- stackApply(species.stack, c(1), sum, na.rm=FALSE)
        threshold <- length(env.names)
        species.layer <- species.layer >= threshold
        return(species.layer)
    } else {
        species.layer <- env.stack[[env.names[1]]]
        col.number <- length(species.list)
        species.stack <- stack()
        for (col.index in 1:col.number) {
            # species.raster <- setValues(species.layer, as.vector(species.list[[col.index]]))
            species.raster <- setValues(species.layer, getValues(species.list[[col.index]]))
            species.stack <- stack(species.stack, species.raster)
        }
        names(species.stack) <- env.names
        return(species.stack)
    }
}

.pickMedian <- function(env.name, env.stack) {
    env.layer <- env.stack[[env.name]]
    env.matrix <- getValues(env.layer)
    env.quantile <- quantile(env.matrix, na.rm=TRUE)
    result.layer <- ((env.layer >= env.quantile[2]) & (env.layer <= env.quantile[4]))
    return(result.layer)
}
