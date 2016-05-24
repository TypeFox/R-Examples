
pfLineartBS<- function(data, particles=1000, plot=FALSE, onlinePlot) {

    # if no data supplied, use default
    if (missing(data)) data <- getPfLineartBSData()

    if (missing(onlinePlot)) {
        useOnline <- FALSE
        onlinePlot <- function() { NULL }
    } else {
        useOnline <- TRUE
        # set up graphics window
        dev.new(width=3,height=3)
        par(mar=c(3,3,1,1),cex=0.8, pch=19, ask=FALSE)
    }

    # more eloquent tests can be added
    stopifnot(nrow(data) > 0,
              ncol(data) == 2,
              colnames(data) == c("x", "y"),
              class(onlinePlot) == "function")

    res <- .Call("pfLineartBS", as.matrix(data),
                 particles,
                 useOnline,
                 onlinePlot,
                 package="RcppSMC")

    if (plot) {
        ## plot 5.1 from vignette / paper
        with(data, plot(x, y, col="red"))
        with(res, lines(Xm, Ym, lty="dashed"))
    }

    invisible(res)
}

# simple convenience function, should probably make the data a
# data component of the package...
getPfLineartBSData <- function() {
    file <- system.file("sampleData", "pf-data.csv", package="RcppSMC")
    dat <- read.table(file, skip=1, header=FALSE, col.names=c("x","y"))
    invisible(dat)
}

pfLineartRange <- function(rrng)
{
   min <- rrng[1]
   max <- rrng[2]

   if(min > 0) {
      rmin = exp(floor(log(min)))
   } else if (min < 0) {
      rmin = -exp(ceiling(log(-min)))
   } else {
      rmin = 0
   }

   if(max > 0) {
      rmax = exp(ceiling(log(max)))
   } else if (max < 0){
      rmax = exp(floor(log(-max)))
   } else {
      rmax = 0;
   }

   invisible(c(rmin,rmax))
}

pfLineartBSOnlinePlot <- function(xm, ym) {
    plot(xm, ym, xlim = pfLineartRange(range(xm)), ylim=pfLineartRange(range(ym)))
    # FIXME sleep time should also be a variable
    Sys.sleep(0.05)
}

