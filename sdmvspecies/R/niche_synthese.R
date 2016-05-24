#!/usr/bin/env Rscript

.gaussianTranslate <- function(factor, range, min, mean) {
    normal.sd <- range/6
    result <- dnorm(factor, mean, normal.sd)

    result <- sqrt(2*pi)*normal.sd*result
    return(result)
}


.linearIncreaseTranslate <- function(factor, range, min, mean) {
    result <- 1/range*(factor-min)
}


.linearDecreaseTranslate <- function(factor, range, min, mean) {
    result <- 1/range*(min+range-factor)
}


.truncatedLinearIncreaseTranslate <- function(factor, range, min, mean) {
    critical.point <- range*2/3+min

    original.factor <- factor
    factor.na <- is.na(original.factor)
    indictor <- original.factor >= critical.point
    indictor[factor.na] <- FALSE

    factor[indictor]<- 1
    indictor <- original.factor < critical.point
    indictor[factor.na] <- FALSE
    new.value <- 1/critical.point*(original.factor-min)
    factor[indictor] <- new.value[indictor]

    return(factor)
}


.truncatedLinearDecreaseTranslate <- function(factor, range, min, mean) {
    critical.point <- range*1/3+min

    original.factor <- factor
    factor.na <- is.na(original.factor)
    indictor <- original.factor <= critical.point
    indictor[factor.na] <- FALSE

    factor[indictor]<- 1
    indictor <- original.factor > critical.point
    indictor[factor.na] <- FALSE
    new.value <- 1/(range-critical.point)*(range+min-original.factor)
    factor[indictor] <- new.value[indictor]

    return(factor)
}


#' nicheSynthese
#'
#' niche synthese method
#'
#' This method mainly implement niche synthese method, for more details see references
#'
#' You can write several paragraphs.
#' @param env.stack a \code{rasterStack} object that contain the environment variable
#' @param config config is a \code{list} or \code{matrix} or \code{data.frame} that contain config info, details see details part
#' @param stack stack is an option that if you want not compose them togethor (result return as a \code{rasterStack}). Default is FALSE
#' @param random.error add random error on cell or not. Default is FALSE
#' @return \code{rasterLayer} or \code{rasterStack} if stack is set to TRUE
#' @references Hirzel, A. H., Helfer, V., & Metral, F. (2001). Assessing habitat-suitability models with a virtual species. Ecological modelling, 145(2), 111-121.
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
#' config <- list(
#'     c("bio1","1",2),
#'     c("bio14", "2", 2),
#'     c("bio5", "3", 1),
#'     c("bio11", "4", 2),
#'     c("bio16", "5", 1)
#' )
#' # call the niche synthsis method
#' species.raster <- nicheSynthese(env.stack, config)
#' # let see the result raster,
#' # you should noticed that it's continue value map not distributin map
#' species.raster
#'
#' # write the map to file, so you can use it latter in GIS software
#' # or further analysis.
#' #
#' #writeRaster(species.raster, "synthese.img", "HFA", overwrite=TRUE)
#'
#' # to make binary distribution map, you should chosee a threshold to make map
#' # see the map then to decide the threshold to binary
#' plot(species.raster)
#' # choice threshold, here we choice 4
#' threshold <- 14
#' # make binary map
#' distribution.map <- species.raster > threshold
#' # plot the map out
#' plot(distribution.map)
nicheSynthese <- function(env.stack, config, stack=FALSE, random.error=FALSE) {
    RESPONSE_METHOD = seq(1,5)

    # check env.stack first
    if (!(class(env.stack) %in% "RasterStack")) {
        stop("env.stack is not a RasterStack object!")
    }
    # check params
    for (config.item in config) {
        #print(param.item)
        if (!(config.item[1] %in% names(env.stack))) {
            stop("params must match with env.stack in layer names")
        }
        if (!(config.item[2] %in% RESPONSE_METHOD)) {
            stop("response method is wrong")
        }
    }
    # TODO:here used mclapply but not given core.number
    species.list <- mclapply(X=config, FUN=.nicheSyntheseMain, env.stack)

    species.matrix <- matrix(unlist(species.list), ncol=length(config), byrow=FALSE)
    if (!stack) {
        species <- apply(species.matrix, 1, sum)
        # calcuate weight
        weight.list <- lapply(X=config, FUN=function(config.item) {return(config.item[3])})
        weight.value <- sum(as.integer(unlist(weight.list)))
        species <- species/weight.value
        # add random error term on condition
        if (random.error) {
            cell.number <- length(species)
            random.error <- runif(cell.number, -0.05, 0.05)
            species <- species + random.error
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
        ncol <- length(config[[1]])
        config.matrix <- matrix(unlist(config), byrow=TRUE, ncol=ncol)
        layer.names <- config.matrix[,1]
        names(species.stack) <- layer.names
        return(species.stack)
    }
}

.nicheSyntheseMain <- function(var, env.stack) {
    predictor.name <- var[1]
    niche.function <- var[2]
    weight <- var[3]
    weight <- as.integer(weight)

    env.layer <- env.stack[[predictor.name]]

    factors.max <- cellStats(env.layer, stat='max', na.rm=TRUE)
    factors.min <- cellStats(env.layer, stat='min', na.rm=TRUE)
    factors.range <- factors.max - factors.min
    factors.mean <- factors.range*0.5 + factors.min

    factor <- getValues(env.layer)

    if (niche.function == "1") {
        result <- .gaussianTranslate(factor, factors.range, factors.min, factors.mean)
    } else if (niche.function == "2") {
        result <- .linearIncreaseTranslate(factor, factors.range, factors.min, factors.mean)
    } else if (niche.function == "3") {
        result <- .linearDecreaseTranslate(factor, factors.range, factors.min, factors.mean)
    } else if (niche.function == "4") {
        result <- .truncatedLinearIncreaseTranslate(factor, factors.range, factors.min, factors.mean)
    } else if (niche.function == "5") {
        result <- .truncatedLinearDecreaseTranslate(factor, factors.range, factors.min, factors.mean)
    }
    result <- result*weight
    return(result)
}
