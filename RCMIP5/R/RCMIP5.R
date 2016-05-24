#' Tools for Manipulating and Summarizing CMIP5 Data
#'
#' Working with CMIP5 data can be tricky, forcing scientists to
#' write custom scripts and programs. The `RCMIP5` package aims
#' to ease this process, providing a standard, robust, and
#' high-performance set of functions to (i) explore what data
#' have been downloaded, (ii) identify missing data, (iii)
#' average (or apply other mathematical operations) across
#' experimental ensembles, (iv) produce both temporal and spatial
#' statistical summaries, and (v) produce
#' easy-to-work-with graphical and data summaries.
#'
#' ...
#'
#' @references Todd-Brown and Bond-Lamberty, 2014: (in prep).
#' @references Taylor et al., 2012:
#'   An overview of CMIP5 and the experiment design, Bulletin of the American
#'   Meteorological Society, 93, 485-498.
#'   \url{http://dx.doi.org/10.1175/BAMS-D-11-00094.1}
#' @import dplyr reshape2 digest abind
#' @docType package
#' @name RCMIP5
NULL

#' The 'cmip5data' class
#'
#' This constructor has two functions. First, given a list, it makes the list
#' a cmip5data-class object (no check is made that the list has appropriate
#' fields though). Second, if given a numeric value(s), it returns sample/
#' example data in the newly constructed object. This is used extensively by
#' the testing code.
#'
#' @param x A list or numeric. If x is a list then the fields are expected to match those of the returned cmip5data object. If x is a numeric sample data is created where the numeric indicates the years of sample data to return.
#' @param lonlat Boolean indicating whether to create lon and lat dimensions
#' @param lonsize Integer size of longitude dimension
#' @param latsize Integer size of latitude dimension
#' @param Z logical. Create Z dimension?
#' @param Zsize integer. Size of Z dimension
#' @param time logical. Create time dimension?
#' @param monthly logical. Monthly (if not, annual) data?
#' @param randomize logical. Random sample data?
#' @param verbose logical. Print info as we go?
#' @return A cmip5data object, which is a list with the following fields:
#'  \item{files}{Array of strings containg the file(s) included in this dataset}
#'  \item{variable}{String containg the variable name described by this dataset}
#'  \item{model}{String containing the model name of this dataset}
#'  \item{experiment}{String containing the experiment name of this dataset}
#'  \item{ensembles}{Array of strings containg the ensemble(s) included in this dataset}
#'  \item{domain}{String containing the domain name of this dataset}
#'  \item{val}{Data frame holding data, with fields lon, lat, Z, time}
#'  \item{valUnit}{String containing the value units}
#'  \item{lon}{Numeric vector containing longitude values; may be \code{NULL}}
#'  \item{lat}{Numeric vector containing latitude values; may be \code{NULL}}
#'  \item{Z}{Numeric vector Z values; may be \code{NULL}}
#'  \item{time}{Numeric vector containing time values; may be \code{NULL}}
#'  \item{dimNames}{Array of strings containing the original (NetCDF) dimension names}
#'  \item{calendarStr}{String defining the calendar type; may be \code{NULL}}
#'  \item{debug}{List with additional data (subject to change)}
#'  \item{provenance}{Data frame with the object's provenance. See \code{\link{addProvenance}}}
#'  \item{numPerYear}{Numeric vector; only present after \code{\link{makeAnnualStat}}}
#'  \item{numYears}{Numeric vector; only present after \code{\link{makeMonthlyStat}}}
#'  \item{numCells}{Numeric vector; only present after \code{\link{makeGlobalStat}}}
#'  \item{filtered}{Logical; only present after \code{\link{filterDimensions}}}
#' @docType class
#' @examples
#' cmip5data(1970)  # produces monthly sample data for year 1970
#' cmip5data(1970:2014)
#' cmip5data(1970:2014, monthly=FALSE)  # annual data
#' cmip5data(1970:2014, randomize=TRUE) # randomized data
#' cmip5data(1970:2014, Z=TRUE)  # four-dimensional data
#' cmip5data(0, time=FALSE)  # sample 'fx' data, two-dimensional
#' cmip5data(list())  # makes this (here empty) list class into 'cmip5data'
#' @export
cmip5data <- function(x=list(),
                      # parameters for making sample data
                      lonlat=TRUE, lonsize=10, latsize=10,
                      Z=FALSE, Zsize=5,
                      time=TRUE, monthly=TRUE,
                      randomize=FALSE, verbose=FALSE) {
    
    # Sanity checks
    stopifnot(is.numeric(c(lonsize, latsize, Zsize)))
    stopifnot(is.logical(c(lonlat, Z, time, monthly, randomize)))
    
    if (is.list(x)) {          # If x is a list then we are done.
        # Just cast it directly to a cmip5data object
        if(verbose) cat("Casting list to cmip5data\n")
        structure(x, class="cmip5data")
        
    } else if(is.numeric(x)) {  # Create sample data
        if(verbose) cat("Creating new cmip5data\n")
        
        # Construct two lists which will be used to create the sample data:
        # ... result and debug.
        
        # result holds the primary data of interest
        result <- list(
            files=NULL,
            variable="var",
            model="model",
            experiment="experiment",
            ensembles="ensemble",
            domain="domain",
            val=NULL,
            valUnit=NULL,
            lon=NA,
            lat=NA,
            Z=NA,
            time=NA,
            dimNames=NULL
        )
        
        debug <- list()
        
        # If this data will have spatial dimensions, construct
        if(lonlat) {
            if(verbose) cat("Adding spatial dimensions\n")
            
            # realistic lon (0 to 360) and lat (-90 to 90) numbers
            result$lon <- 360/lonsize * c(0:(lonsize-1))  + 360/lonsize/2
            result$lat <- 180/latsize * c(0:(latsize-1)) - 90 + 180/latsize/2
            result$dimNames=c("lon", "lat")
            debug$lonUnit <- "degrees_east"
            debug$latUnit <- "degrees_north"
        } else {
            result$dimNames <- c(NA, NA)
        }
        
        # If this data will have Z dimension, construct
        if(Z) {
            if(verbose) cat("Adding Z dimensions\n")
            
            result$Z <- c(0:(Zsize-1))
            result$dimNames <- c(result$dimNames, "Z")
            debug$ZUnit <- "m"
        } else {
            result$dimNames <- c(result$dimNames, NA)
        }
        
        # If this data will have time dimension, construct
        if(time) {
            if(verbose) cat("Adding time dimensions\n")
            
            years <- x
            ppy <- ifelse(monthly, 12, 1)  # periods per year
            result$calendarStr <- "360_day"
            debug$timeFreqStr <- ifelse(monthly, "mon", "yr")
            debug$startYr <- years[1]
            debug$calendarStr <- "360_day"
            debug$timeUnit <- paste0("days since ",years[1],"-01-01")
            
            if(monthly) {
                # '+15' initalizes all time stamps to be middle of the month
                debug$timeRaw <- (360/ppy*c(0:(length(years)*ppy-1) )+15)    
                result$time <- debug$timeRaw/360+min(years)    
            } else {
                debug$timeRaw <- result$time <- years
            }
            
            # convert day based calandar to year based
            result$dimNames <- c(result$dimNames, "time")
        } else { # no time
            result$dimNames <- c(result$dimNames, NA)
            result$domain <- "fx"
        }
        
        # Make data frame, fill it with fake data, wrap as tbl_df
        result$val <- expand.grid(lon=result$lon, lat=result$lat,
                                  Z=result$Z, time=result$time)
        if(randomize) {
            result$val$value <- runif(n=nrow(result$val))
        } else {
            result$val$value <- 1
        }      
        result$val <- tbl_df(result$val)
        
        result$valUnit <- "unit"
        result$debug <- debug
        
        # Change any NA dimension (was needed for expand.grid above) to NULL
        if(all(is.na(result$lon))) result$lon <- NULL
        if(all(is.na(result$lat))) result$lat <- NULL
        if(all(is.na(result$Z))) result$Z <- NULL
        if(all(is.na(result$time))) result$time <- NULL
        
        # Add debug info and set class
        result <- structure(result, class="cmip5data")
        
        # Initialize provenance and return
        addProvenance(result, "Dummy data created")
    } else {
        stop("Don't know what to do with this class of parameter")
    }
}

#' Print a 'cmip5data' class object.
#'
#' @param x A \code{\link{cmip5data}} object
#' @param ... Other parameters passed to cat
#' @details Prints a one-line summary of the object
#' @method print cmip5data
#' @export
#' @keywords internal
print.cmip5data <- function(x, ...) {
    
    if(is.null(x$variable)) {
        cat("(Empty cmip5data object)")
        return()
    }
    
    ansStr <- paste0('CMIP5: ', x$variable, ", ", x$model, " ", x$experiment)
    
    spaceStr <- paste(length(x$lon), "x", length(x$lat), "x", length(x$Z))
    ansStr <- paste0(ansStr, ", ", spaceStr)
    
    timeStr <- "no time"
    if(!is.null(x$time) & length(x$time) > 0) {
        timeStr <- paste(floor(min(x$time, na.rm=TRUE)), "to",
                         floor(max(x$time, na.rm=TRUE)))
    }
    ansStr <- paste0(ansStr, ", ", timeStr)
    
    if(!is.null(x$ensembles)) {
        ansStr <- paste0(ansStr, ", from ", length(x$ensembles), " ",
                         ifelse(length(x$ensembles)==1, "ensemble", "ensembles"))
    }
    
    cat(ansStr, "\n")
    cat(...)
} # print.cmip5data

#' Summarize a 'cmip5data' class object.
#'
#' @param object A \code{\link{cmip5data}} object
#' @param ... ignored
#' @details Prints a short summary of the object.
#' @return A summary structure of the object.
#' @method summary cmip5data
#' @export
#' @keywords internal
summary.cmip5data <- function(object, ...) {
    
    ans <- list()
    class(ans) <- "summary.cmip5data"
    
    # cmip5 objects should always have the following defined:
    ans$variable <- object$variable
    ans$valUnit <- object$valUnit
    ans$domain <- object$domain
    ans$model <- object$model
    ans$experiment <- object$experiment
    ans$ensembles <- object$ensembles
    
    ans$type <- "CMIP5 data"
    
    #    if(grepl('makeAnnualStat', rev(object$provenance$caller)[1])) {
    #        ans$type <-  paste(ans$type, "annual summary [",
    #                           paste(unique(object$numPerYear), collapse=', '),
    #                           "]") 
    #    }
    
    if(!is.null(object$numPerYear)) {
        ans$type <-  paste0(ans$type, " (annual summary of ", mean(object$numPerYear), " times)")  
    }
    if(!is.null(object$numYears)) {
        ans$type <- paste0(ans$type, " (monthly summary of ", mean(object$numYears), " years)")
    } 
    if(!is.null(object$numCells)) {
        ans$type <- paste0(ans$type, " (spatial summary of ", object$numCells, " cells)")
    } 
    if(!is.null(object$numZs)) {
        ans$type <- paste0(ans$type, " (Z summary of ", object$numZs, " levels)")
    } 
    
    if(!is.null(object$filtered)) {
        ans$type <- paste(ans$type, "(filtered)")
    }
    
    #    if(!is.null(object$area)) {
    #        ans$type <- paste(ans$type, "(regridded)")
    #    }
    
    ans$spatial <- paste0("lon [", length(object$lon),
                          "] lat [", length(object$lat),
                          "] Z [", length(object$Z), "]")
    
    ans$time <- paste0(object$debug$timeFreqStr, " [", length(object$time), "] ", object$debug$timeUnit)
    ans$size <- as.numeric(object.size(object))
    ans$valsummary <- c(min(object$val$value, na.rm=TRUE),
                        mean(object$val$value, na.rm=TRUE),
                        max(object$val$value, na.rm=TRUE))
    ans$provenance <- object$provenance
    
    return(ans)
} # summary.cmip5data

#' Print the summary for a 'cmip5data' class object.
#'
#' @param x A \code{\link{cmip5data}} object
#' @param ... Other parameters passed to cat
#' @details Prints a one-line summary of the object
#' @method print summary.cmip5data
#' @export
#' @keywords internal
print.summary.cmip5data <- function(x, ...) {
    cat(x$type, "\n")
    cat("Variable: ", x$variable, " (", x$valUnit, ") from model ", x$model, "\n", sep="")
    cat("Data range: ", round(x$valsummary[1], 2), "-", round(x$valsummary[3], 2),
        "  Mean: ", round(x$valsummary[2], 2), "\n", sep="")
    cat("Experiment:", x$experiment, "-", length(x$ensembles), "ensemble(s)\n")
    cat("Spatial dimensions:", x$spatial, "\n")
    cat("Time dimension:", x$time, "\n")
    cat("Size:", format(round(x$size/1024/1024, 1), nsmall=1), "MB\n")
    cat("Provenance has", nrow(x$provenance), "entries\n")
} # print.summary.cmip5data

#' Convert a cmip5data object to a data frame
#'
#' @param x A \code{\link{cmip5data}} object
#' @param ... Other parameters
#' @param originalNames logical. Use original dimension names from file?
#' @return The object converted to a data frame
#' @export
#' @keywords internal
as.data.frame.cmip5data <- function(x, ..., originalNames=FALSE) {
    # Suppress stupid NOTEs from R CMD CHECK
    lon <- lat <- Z <- time <- NULL
    dplyr::arrange(x$val, lon, lat, Z, time)
} # as.data.frame.cmip5data

#' Convert a cmip5data object to an array
#'
#' @param x A \code{\link{cmip5data}} object
#' @param ... Other parameters
#' @param drop logical. Drop degenerate dimensions?
#' @return The object converted to an array
#' @export
#' @keywords internal
as.array.cmip5data <- function(x, ..., drop=TRUE) {
    
    dimList <- c(length(unique(x$val$lon)),
                 length(unique(x$val$lat)),
                 length(unique(x$val$Z)),
                 length(unique(x$val$time)))
    
    # Remove degenerate dimensions
    if(drop) {
        dimList <- dimList[!dimList %in% 1]        
    }
    
    # Suppress stupid NOTEs from R CMD CHECK
    lon <- lat <- Z <- time <- NULL
    
    # Note we sort data frame before converting to array!
    array(dplyr::arrange(x$val, lon, lat, Z, time)$value, dim=dimList)
} # as.array.cmip5data

#' Make package datasets and write them to disk.
#'
#' @param path root of directory tree
#' @param maxSize max size (in MB) of dataset to write
#' @param outpath directory to write to
#' @details Writes all available ensembles to disk as Rdata files, subject to
#' a maximum size parameter (CRAN says keep sample data < 5MB).
#' @note This is an internal RCMIP5 function and not exported.
#' @keywords internal
makePackageData <- function(path="./sampledata", maxSize=Inf, outpath="./data") {
    if(!file.exists(outpath)) dir.create(outpath)
    stopifnot(file.exists(outpath))
    datasets <- getFileInfo(path)
    if(is.null(datasets)) return()
    
    for(i in 1:nrow(datasets)) {
        cat("-----------------------\n", datasets[i, "filename"], "\n")
        d <- with(datasets[i,],
                  loadEnsemble(variable, model, experiment, ensemble, path=path, verbose=T)
        )
        print(object.size(d), units="MB")
        if(object.size(d)/1024/1024 <= maxSize) {
            objname <- gsub("_[0-9]{4,}-[0-9]{4,}.nc$", "", basename(d$files[1])) # strip dates
            assign(objname, d)
            cat("Writing", objname, "\n")
            save(list=objname, file=paste0(outpath, "/", objname, ".rda"))
        } else {
            cat("Too big; skipping\n")
        }
    }
} # makePackageData
