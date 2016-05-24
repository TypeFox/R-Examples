###########################################################################
##                                                                       ##
## plot.cma      - plots close modern analogue results                   ##
##                                                                       ##
## Created       : 27-May-2006                                           ##
## Author        : Gavin Simpson                                         ##
## Version       : 0.1                                                   ##
## Last modified : 27-May-2006                                           ##
##                                                                       ##
## ARGUMENTS:                                                            ##
## x             - object of class cma.                                  ##
## method        - the method to be used to separate coincident points.  ##
##                 default method '"overplot"' causes such points to be  ##
##                 overplotted, but it is also possible to specify       ##
##                 "jitter" to jitter the points, or '"stack"' have      ##
##                 coincident points stacked.  The last method only      ##
##                 makes sense for very granular data.character. As per  ##
##                 arguments for stripchart.                             ##
## jitter        - when 'method="jitter"' is used, 'jitter' gives the    ##
##                 amount of jittering applied.                          ##
## vertical      - when vertical is 'TRUE' the plots are drawn           ##
##                 vertically rather than the default horizontal.        ##
## main, ylab, xlab,                                                     ##
##               - graphical parameters.                                 ##
## ...           - additional arguments passed to stripchart             ##
##                                                                       ##
###########################################################################
plot.cma <- function(x, method = c("overplot", "jitter", "stack"),
                     jitter = 0.1, vertical = FALSE, draw.quant = TRUE,
                     xlab = NULL, ylab = "", main = "", cex.axis = NULL,
                     ..., col.quant = "red", lty.quant = "dashed")
  {
    if (!inherits(x, "cma")) 
      stop("use only with \"cma\" objects")
    if(is.null(cex.axis))
      opar <- par(mar = c(5, 5, 4, 2) + 0.1, las = 1)
    else
      opar <- par(mar = c(5, 5, 4, 2) + 0.1, las = 1,
                  cex.axis = cex.axis)
    on.exit(par(opar))
    x <- summary(x)
    dat <- as.vector(x$distances)
    dims <- dim(x$distances)
    NAs <- !is.na(dat)
    dat <- dat[NAs]
    method <- match.arg(method)
    groups <- gl(dims[2], dims[1], labels = colnames(x$distances))
    groups <- factor(groups[NAs], exclude = TRUE)
    if(is.null(xlab))
      xlab <- paste("Dissimilarity <", round(x$cutoff, 4))
    stripchart(dat ~ groups, method = method, vertical = vertical,
               jitter = jitter, main = main, xlab = xlab, ylab = ylab,
               ...)
    if(draw.quant) {
      sel <- x$quant <= x$cutoff
      if(any(sel)) {
        quant <- x$quant[sel]
        abline(v = quant, col = col.quant, lty = lty.quant)
        suffix <- rep("th", times = length(x$prob[sel]))
        suffix[which(x$prob == 0.01)] <- "st"
        quant.title <- paste(100 * x$prob[sel], suffix, sep = "")
        axis(side = 3, at = quant, labels = quant.title)
      } else {
        warning("No quantiles within 'x$cutoff'. 'draw.quant' ignored.")
      }
    }
    invisible(list(distances = dat, groups = groups))
  }
