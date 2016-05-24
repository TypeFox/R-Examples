#' Compute annual statistic of a variable
#'
#' Most CMIP5 data are monthly, and we frequently want to summarize these to annual
#' numbers. This function does that (although annual files also occur, and will be
#' handled as well). The default statistic is \link{mean}, but any summary
#' function that returns a numeric result can be used.
#'
#' @param x A \code{\link{cmip5data}} object
#' @param verbose logical. Print info as we go?
#' @param FUN function. Function to apply across months of year
#' @param ... Other arguments passed on to \code{FUN}
#' @return A \code{\link{cmip5data}} object, whose \code{val} field is the annual
#' mean of the variable. A \code{numMonths} field is also added
#' recording the number of months averaged for each year.
#' @details The stat function is calculated for all combinations of lon,
#' lat, and Z (if present).
#' @examples
#' d <- cmip5data(1970:1975)   # sample data
#' makeAnnualStat(d)
#' summary(makeAnnualStat(d))
#' summary(makeAnnualStat(d, FUN=sd))
#' @seealso \code{\link{makeZStat}} \code{\link{makeGlobalStat}} \code{\link{makeMonthlyStat}}
#' @export
makeAnnualStat <- function(x, verbose=FALSE, FUN=mean, ...) {
    
    # Sanity checks
    stopifnot(class(x)=="cmip5data")
    
    stopifnot(length(verbose)==1 & is.logical(verbose))
    stopifnot(length(FUN)==1 & is.function(FUN))
    
    # Main computation code
    timer <- system.time({  # time the main computation, below
        # Put data in consistent order BEFORE overwriting time
        x$val <- group_by(x$val, lon, lat, Z, time) %>%
            arrange()
        
        x$val$year <- floor(x$val$time)   
        
        # Suppress stupid NOTEs from R CMD CHECK
        lon <- lat <- Z <- time <- year <- value <- `.` <- NULL
        
        # Put data in consistent order and compute
        
        # Instead of "summarise(value=FUN(value, ...))", we use the do()
        # call below, because the former doesn't work (as of dplyr 0.3.0.9000):
        # the ellipses cause big problems. This solution thanks to Dennis
        # Murphy on the manipulatr listesrv.
        x$val <- group_by(x$val, lon, lat, Z, year) %>%
            do(data.frame(value = FUN(.$value, ...))) %>%
            ungroup()
        x$val$time <- x$val$year
        x$val$year <- NULL
        
        # dplyr doesn't (yet) have a 'drop=FALSE' option, and the summarise
        # command above may have removed some lon/lat combinations
        if(length(unique(x$val$lon)) < length(x$lon) |
               length(unique(x$val$lat)) < length(x$lat)) {
            if(verbose) cat("Replacing missing lon/lat combinations\n")
            
            # Fix this by generating all lon/lat pairs and combining with answer
            full_data <- tbl_df(expand.grid(lon=x$lon, lat=x$lat))
            x$val <- left_join(full_data, x$val, by=c("lon", "lat"))
        }
    }) # system.time
    
    if(verbose) cat('\nTook', timer[3], 's\n')
    
    x$numPerYear <- as.data.frame(table(floor(x$time)))$Freq
    x$time <- unique(floor(x$time))
    x$debug$timeFreqStr <- "years (summarized)"
    addProvenance(x, paste("Calculated", 
                           paste(deparse(substitute(FUN)), collapse="; "),
                           "for years", min(x$time), "-", max(x$time)))
} # makeAnnualStat
