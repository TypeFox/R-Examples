#-----------------------------------------------------------------------------#
#                                                                             #
#                     QUALITY CONTROL STATISTICS IN R                         #
#                                                                             #
#  An R package for statistical in-line quality control.                      #
#                                                                             #
#  Written by: Miguel A. Flores Sánchez                                       #
#              Student Master of Statistical Techniques                       #
#              University of The Coruña, SPAIN                                #
#              mflores@outlook.com                                            #
#                                                                             #
#-----------------------------------------------------------------------------#
#-------------------------------------------------------------------------
# cusum chart
#-------------------------------------------------------------------------
##' Function to plot cusum chart
##'
##' This function is used to compute statistics required by the cusum chart.
##'
##' @param x   an R object (used to select the method). See details.
##' @param ... arguments passed to or from methods.
##' @export
##' @examples
##' 
##' library(qcr)
##' data(pistonrings)
##' attach(pistonrings)
##' res.qcd <- qcd(pistonrings, type.data = "dependence")
##' res.qcs <- qcs.cusum(res.qcd, type = "cusum")
##' summary(res.qcs)
##' plot(res.qcs)
##' 
qcs.cusum <- function(x, ...) {
  UseMethod("qcs.cusum")
}

##' @rdname qcs.cusum
##' @method qcs.cusum default
##' @inheritParams qcd
##' @param sizes a value or a vector of values specifying the sample sizes
##' associated with each group.
##' @param center a value specifying the center of group statistics or the
##' ''target'' value of the process.
##' @param std.dev a value or an available method specifying the within-group
##' standard deviation(s) of the process. \cr Several methods are available for
##' estimating the standard deviation.
##' @param decision.interval A numeric value specifying the number of standard
##' errors of the summary statistics at which the cumulative sum is out of
##' control.
##' @param se.shift The amount of shift to detect in the process, measured in
##' standard errors of the summary statistics.
##' @param plot a logical value indicating should be plotted.
##' @export
##' 
qcs.cusum.default <- function(x, var.index  =  1, sample.index  =  2,
                          covar.index  =  NULL, covar.names  =  NULL,
                          data.name = NULL,
                          sizes = NULL,
                          center = NULL, std.dev = NULL,
                          decision.interval =  5, 
                          se.shift = 1, plot = FALSE, ...)
{
  obj<-qcd(data = x, var.index = var.index, sample.index = sample.index,
           covar.index = covar.index, covar.names = covar.names,
           data.name = data.name, sizes = sizes, type.data = "dependence")
  
  result<-qcs.cusum.qcd(x = obj, center = center, std.dev = std.dev, 
                       decision.interval = decision.interval, 
                       se.shift = se.shift, plot = plot)
  
  return(result)
}


##' @rdname  qcs.cusum
##' @method qcs.cusum qcd
##' @inheritParams qcs.cusum.default
##' @export
##' 
qcs.cusum.qcd <- function(x, center = NULL,
                         std.dev = NULL, 
                         decision.interval  =  5, se.shift = 1, plot = FALSE, ...) {
  #.........................................................................
  if(is.null(x) || !inherits(x, "qcd"))
    stop("data must be an objects of class (or extending) 'qcd'")
  sizes <- x$sizes
  type.data <- "dependence" 

  std <- if(any(sizes==1)) "xbar.one" else "xbar"
  if(is.null(std.dev)) 
  { std.dev <- switch(std, 
                      "xbar" = { if(any(sizes > 25)) "RMSDF"
                                 else                "UWAVE-R" },
                      "xbar.one" = "MR")
  }  
  
  qcs<-qcs(x = x$x, sample.index = x$sample, sizes = sizes,
           center = center, std.dev = std.dev, type = "cusum",
           decision.interval = decision.interval, se.shift = se.shift, type.data = type.data)
  
  center <- qcs$center
  cusum <- qcs$statistics
  std.dev <- qcs$std.dev
  sizes <- qcs$sizes
  limits <- qcs$limits
  violations <- qcs$violations
  pos <- qcs$pos
  neg  <- qcs$neg
  decision.interval <- qcs$decision.interval
  se.shift <- qcs$se.shift
  
  statistics <- data.frame(cusum)
  m <- length(x)
  sample <- x$sample
  
  if (m > 3) {
    new.x <- x[, -c(1, 2, length(x))]
    cov <- apply(new.x, 2, function(x) unlist(lapply(split(x, sample), unique)))
    statistics <- data.frame(cusum, cov)
  }
  
  row.names(statistics) <- unique(x$sample)
  data.name <- attr(x, "data.name")
  result <- list(qcd  =  x, type  =  "cusum", statistics  =  statistics,
                 center  =  center, std.dev  =  std.dev,
                 limits  =  limits, 
                 sizes  =  sizes, data.name  =  data.name,
                 violations  =  violations, pos = pos, neg = neg,
                 decision.interval = decision.interval, 
                 se.shift = se.shift)
  
  oldClass(result) <- c("qcs.cusum", "qcs")
  
  if(plot) plot(result, ...)
  
  return(result)
  #.........................................................................
} # qcs.cusum.qcd