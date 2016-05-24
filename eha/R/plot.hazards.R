plot.hazards <- function(x,
                         fn = c("cum", "surv", "log", "loglog"),
                         fig = TRUE,
                         xlim = NULL,
                         ylim = NULL,
                         main = NULL,
                         xlab = NULL,
                         ylab = NULL,
                         ...
                         ){
    ## x is of type 'hazards', which is a data.frame with two or three components,
    ## 'hazard' (e.g. cumulative...), 'time', and optionally 'strata'
    ## as components. The first column contains
    ## risktimes, and the second column the corresponding 'hazard atoms'.

            ##cat("length(x$time) = ", length(x$time), "\n") 
            ##cat("length(x$hazard) = ", length(x$hazard), "\n") 
    fn <- fn[1]
    
    if (is.null(x)){
        cat("Must be fixed in plot.hazdata!\n")
        return(NULL)
    }
        
    if (is.null(x$hazard) || is.null(x$time))
        stop("'x' is of wrong type.")

    if (!(fn %in% c("cum", "surv", "log", "loglog")))
        stop(paste(fn, "is an illegal value of 'fn'"))

    if (!is.null(x$strata)){
        ns <- length(levels(x$strata))
    }else{
        ns <- 1
    }
    max.x <- max(x$time)
    min.x <- min(x$time)
    max.y <- max(x$hazard)
    min.y <- min(x$hazard)
    
    if (fn == "loglog"){
        x$time <- log(x$time)
        max.x <- log(max.x)
        min.x <- log(min.x)
        x$hazard <- log(x$hazard)
        max.y <- log(max.y)
        min.y <- log(min.y)
    }else if (fn == "log"){
        x$hazard <- log(x$hazard)
        max.y <- log(max.y)
        min.y <- log(min.y)
    }else if (fn == "surv"){
        x$hazard <- exp(-x$hazard)
        y.max <- 1
        y.min <- 0
    }

    if (is.null(xlim)) {
        xlim <- c(min.x, max.x)
        ##if (fn != "loglog") xlim[1] <- 0.0
    }
    if (length(xlim) != 2) stop("'xlim' must be a vector of length two")

    if (is.null(ylim)){
        ylim <- c(min.y, max.y)
        if ( (fn == "surv") || (fn == "cum") ) ylim[1] <- 0.0
    }
    if (length(ylim) != 2) stop("'ylim' must be a vector of length two")

    if (fig){
        if (ns == 1){
            ##cat("length(x$time) = ", length(x$time), "\n") 
            ##cat("length(x$hazard) = ", length(x$hazard), "\n") 
            plot(x$time, x$hazard, type = "s",
             xlim = xlim, ylim = ylim,
             xlab = xlab, ylab = ylab, main = main, lty = 1, ...)
        }else{
            stra <- as.integer(x$strata)
            plot(x$time[stra == 1], x$hazard[stra == 1], type = "s",
                 xlim = xlim, ylim = ylim,
                 xlab = xlab, ylab = ylab, main = main, lty = 1, ...)
                 
            
            for (i in 2:ns){
                lines(x$time[stra == i], x$hazard[stra == i], type = "s", lty = i)
            }
        }
        abline(h = 0)
    }
    invisible(list(x = x, fn = fn))
}

