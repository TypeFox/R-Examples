#########################################################
#### AUTHOR:     Arnost Komarek                      ####
####             (2004)                              ####
####                                                 ####
#### FILE:       densplot2.R                         ####
####                                                 ####
#### FUNCTIONS:  densplot2                           ####
#########################################################

### ======================================
### densplot2
### ======================================
## Slightly modified 'densplot' function of 'coda' library
## 05/05/2004
##
## argument 'plot' added to indicate whether a plot is to be created directly
## or a list with data.frames for future plotting is to be returned
## and some other plotting arguments allowed to be changed by a user
##
densplot2 <-
function (x, plot = TRUE, show.obs = FALSE, bwf, bty = "n", main = "", xlim, ylim, xlab, ylab, ...) 
{
    xx <- as.matrix(x)
    toplot <- list()
    for (i in 1:nvar(x)) {
        y <- xx[, i, drop = TRUE]
        if (missing(bwf)) 
            bwf <- function(x) {
                x <- x[!is.na(as.vector(x))]
                return(1.06 * min(sd(x), IQR(x)/1.34) * length(x)^-0.2)
            }
        bw <- bwf(y)
        width <- 4 * bw
        if (max(abs(y - floor(y))) == 0 || bw == 0) 
            hist(y, prob = TRUE, main = main, ...)
        else {
            scale <- "open"
            if (max(y) <= 1 && 1 - max(y) < 2 * bw) {
                if (min(y) >= 0 && min(y) < 2 * bw) {
                  scale <- "proportion"
                  y <- c(y, -y, 2 - y)
                }
            }
            else if (min(y) >= 0 && min(y) < 2 * bw) {
                scale <- "positive"
                y <- c(y, -y)
            }
            else scale <- "open"
            dens <- density(y, width = width)
            if (scale == "proportion") {
                dens$y <- 3 * dens$y[dens$x >= 0 & dens$x <= 
                  1]
                dens$x <- dens$x[dens$x >= 0 & dens$x <= 1]
            }
            else if (scale == "positive") {
                dens$y <- 2 * dens$y[dens$x >= 0]
                dens$x <- dens$x[dens$x >= 0]
            }
            if (missing(ylim)) ylim <- c(0, max(dens$y))
            if (missing(xlab)) xlab <- paste("N =", niter(x), "  Bandwidth =", formatC(dens$bw))
            if (missing(ylab)) ylab <- ""
            if (missing(xlim)) xlim <- NULL
            if (plot){
                plot(dens, main = main, type = "l", bty = bty, xlab = xlab, ylab = ylab, xlim = xlim, ylim = ylim, ...)
               if (show.obs) 
                   lines(y[1:niter(x)], rep(max(dens$y)/100, niter(x)), 
                     type = "h")
            }                
            else
               toplot[[i]] <- data.frame(x = dens$x, y = dens$y)
        }
        if (!is.null(varnames(x)) && is.null(list(...)$main)) 
            title(paste("Density of", varnames(x)[i]))
    }
    if (plot) return(invisible(x))
    else      return(toplot)
}
