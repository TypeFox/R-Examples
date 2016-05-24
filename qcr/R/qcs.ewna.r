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
# ewma chart
#-------------------------------------------------------------------------
##' Function to plot ewma chart
##'
##' This function is used to compute statistics required by the ewma chart.
##'
##' @param x   an R object (used to select the method). See details.
##' @param ... arguments passed to or from methods.
##' @export

qcs.ewma <- function(x, ...) {
  UseMethod("qcs.ewma")
}

##' Function to plot ewma chart
##'
##' This function is used to compute statistics required by the ewma chart.
##'
##' @rdname qcs.ewma
##' @method qcs.ewma default
##' @inheritParams qcd
##' @param center a value specifying the center of group statistics or the
##' ''target'' value of the process.
##' @param std.dev  a value or an available method specifying the within-group standard
##' deviation(s) of the process. Several methods are available for estimating the
##' standard deviation in case of a continuous process variable.
##' @param nsigma  a numeric value used to compute control limits, specifying the
##' number of standard deviations.
##' @param lambda the smoothing parameter \eqn{0 \le \lambda \le 1}{0 <= lambda
##' <= 1}
##' @param plot a logical value indicating should be plotted.
##' @export
qcs.ewma.default <- function(x, var.index  =  1, sample.index  =  2,
                          covar.index  =  NULL, covar.names  =  NULL,
                          data.name = NULL,
                          sizes = NULL,
                          center = NULL, std.dev = NULL,
                          nsigma  =  3, lambda=0.2, plot = FALSE, ...)
{
  obj<-qcd(data = x, var.index = var.index, sample.index = sample.index,
           covar.index = covar.index, covar.names = covar.names,
           data.name = data.name, sizes = sizes, type.data = "dependence")
  
  result<-qcs.ewma.qcd(x = obj, center = center, std.dev = std.dev, 
                       nsigma = nsigma, lambda = lambda, plot = plot)
  
  return(result)
}


##' @rdname  qcs.ewma
##' @method qcs.ewma qcd
##' @inheritParams qcs.ewma.default
##' @export
qcs.ewma.qcd <- function(x, center = NULL,
                         std.dev = NULL, 
                         nsigma  =  3, lambda = 0.2, plot = FALSE, ...) {
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
           center = center, std.dev = std.dev, type = "ewma",
           conf.nsigma = nsigma, lambda = lambda, type.data = type.data)
  
  center <- qcs$center
  ewma <- qcs$statistics
  std.dev <- qcs$std.dev
  sizes <- qcs$sizes
  limits <- qcs$limits
  violations <- qcs$violations
  xx <- qcs$x
  y  <- qcs$y
  sigma <- qcs$sigma
  lambda <- qcs$lambda
  nsigma  <-  qcs$nsigma
  
  statistics <- data.frame(ewma)
  m <- length(x)
  sample <- x$sample
  
  if (m > 3) {
    new.x <- x[, -c(1, 2, length(x))]
    cov <- apply(new.x, 2, function(x) unlist(lapply(split(x, sample), unique)))
    statistics <- data.frame(ewma, cov)
  }
  
  row.names(statistics) <- unique(x$sample)
  data.name <- attr(x, "data.name")
  result <- list(qcd  =  x, type  =  "ewma", statistics  =  statistics,
                 center  =  center, std.dev  =  std.dev,
                 limits  =  limits, nsigma  =  nsigma,
                 sizes  =  sizes, data.name  =  data.name,
                 violations  =  violations, x = xx, y = y, lambda = lambda, 
                 sigma = sigma)
  
  oldClass(result) <- c("qcs.ewma", "qcs")
  
  if(plot) plot(result, ...)
  
  return(result)
  #.........................................................................
} # qcs.ewma.qcd
