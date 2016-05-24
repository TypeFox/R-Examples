#' Compute monthly statistic of a variable
#' 
#' We frequently want to summarize CMIP5 data by month, e.g. to understand how
#' air temperature varies over the year for a particular data range. This function 
#' does that for monthly data. The default statistic is \link{mean}, but any 
#' summary function that returns a numeric result can be used.
#'
#' @param x A \code{\link{cmip5data}} object
#' @param verbose logical. Print info as we go?
#' @param FUN function. Function to apply across months of year
#' @param ... Other arguments passed on to \code{FUN}
#' @return A \code{\link{cmip5data}} object, whose \code{val} field is the monthly
#' mean of the variable. A \code{numYears} field is also added
#' recording the number of years averaged for each month.
#' @details The stat function is calculated for all combinations of lon,
#' lat, and Z (if present).
#' @seealso \code{\link{makeAnnualStat}} \code{\link{makeZStat}} \code{\link{makeGlobalStat}}
#' @examples
#' d <- cmip5data(1970:1975)   # sample data
#' makeMonthlyStat(d)
#' summary(makeMonthlyStat(d))
#' summary(makeMonthlyStat(d, FUN=sd))
#' @export
makeMonthlyStat <- function(x, verbose=FALSE, FUN=mean, ...) {
    
    # Sanity checks
    stopifnot(class(x)=="cmip5data")
    stopifnot(is.null(x$numYears))
    stopifnot(x$debug$timeFreqStr=="mon")
    stopifnot(length(verbose)==1 & is.logical(verbose))
    stopifnot(length(FUN)==1 & is.function(FUN))
    
    # Main computation code
    timer <- system.time({ # time the main computation, below
        # Put data in consistent order BEFORE overwriting time
        x$val <- group_by(x$val, lon, lat, Z, time) %>%
            arrange()
        
        monthIndex <- floor((x$val$time %% 1) * 12 + 1)
        x$val$month <- monthIndex  
        
        # Suppress stupid NOTEs from R CMD CHECK
        lon <- lat <- Z <- time <- month <- value <- `.` <- NULL
        
        # Put data in consistent order and compute
        
        # Instead of "summarise(value=FUN(value, ...))", we use the do()
        # call below, because the former doesn't work (as of dplyr 0.3.0.9000):
        # the ellipses cause big problems. This solution thanks to Dennis
        # Murphy on the manipulatr listesrv.
        x$val <- group_by(x$val, lon, lat, Z, month) %>%
            do(data.frame(value = FUN(.$value, ...))) %>%
            ungroup()
        x$val$time <- x$val$month
        x$val$month <- NULL
        
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
    
    # Finish up
    x$numYears <- as.data.frame(table(floor(monthIndex)))$Freq
    x$timeUnit <- "months (summarized)"
    x$time <- 1:12
    addProvenance(x, paste("Calculated", 
                           paste(deparse(substitute(FUN)), collapse="; "),
                           "for months", min(x$time), "-", max(x$time)))
} # makeMonthlyStat
