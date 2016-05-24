#-----------------------------------------------------------------------------#
#                                                                             #
#                     QUALITY CONTROL STATISTICS IN R                         #
#                                                                             #
#  An R package for statistical in-line quality control.                      #
#                                                                             #
#  Written by: Miguel A. Flores S?nchez                                       #
#              Student Master of Statistical Techniques                       #
#              University of The Coru?a, SPAIN                                #
#              mflores@outlook.com                                            #
#                                                                             #
#-----------------------------------------------------------------------------#
#-------------------------------------------------------------------------
# xbar chart
#-------------------------------------------------------------------------
##' Function to plot Shewhart xbar chart
##'
##' This function is used to compute statistics required by the xbar chart.
##'
##' @param x   an R object (used to select the method). See details.
##' @param ... arguments passed to or from methods.
##' @export
##' @references Montgomery, D.C. (2000)
##' @examples
##' 
##' ##
##' ##  Continuous data 
##' ##
##'library(qcr)
##'data(pistonrings)
##'str(pistonrings)
##'pistonrings.qcd<-qcd(pistonrings)
##'
##'class(pistonrings.qcd)
##'
##'res.qcs <- qcs.xbar(pistonrings.qcd)
##'plot(res.qcs,title="Control Chart Xbar for pistonrings I")
##'summary(res.qcs)  
##'
##'res.qcd <- state.control(res.qcs)
##'res.qcs <- qcs.xbar(res.qcd)
##'plot(res.qcs,title="Control Chart Xbar for pistonrings II")
##'summary(res.qcs)  
##'
##'res.qcd <- state.control(res.qcs)
##'res.qcs <- qcs.xbar(res.qcd)
##'plot(res.qcs,title="Control Chart Xbar for pistonrings III")
##'summary(res.qcs)  
##'
##'x <- droplevels(pistonrings[1:125,])
##'y <- droplevels(pistonrings[126:200,])
##'
##'res.qcs <- qcs.xbar(x, data.name="Control Chart Xbar for pistonrings")
##'plot(res.qcs)
##'
##'res.qcs <- qcs.add(x = res.qcs, value = y[,c(1,2)]) 
##'plot(res.qcs)
##'summary(res.qcs)
##'
##'
##'res.qcs <- qcs.xbar(pistonrings.qcd, std.dev="UWAVE-SD")
##'class(res.qcs)
##'plot(res.qcs,title="Control Chart Xbar for pistonrings (UWAVE-SD)")
##'summary(res.qcs)  
##'

qcs.xbar <- function(x, ...) {
  UseMethod("qcs.xbar")
}

##' @rdname qcs.xbar
##' @method qcs.xbar default
##' @inheritParams qcd
##' @param center a value specifying the center of group statistics or the
##' ''target'' value of the process.
##' @param std.dev  a value or an available method specifying the within-group standard
##' deviation(s) of the process. Several methods are available for estimating the
##' standard deviation in case of a continuous process variable.
##' @param conf.nsigma  a numeric value used to compute control limits, specifying the
##' number of standard deviations (if \code{conf.nsigma} > 1) or the confidence level (if 0
##' < \code{conf.nsigma} < 1).
##' @param limits a two-values vector specifying control limits.
##' @param plot a logical value indicating should be plotted.
##' @export
##' 
qcs.xbar.default <- function(x, var.index  =  1, sample.index  =  2,
                             covar.index  =  NULL, covar.names  =  NULL,
                             data.name = NULL,
                             sizes  =  NULL,
                             center = NULL,
                             std.dev  =  c("UWAVE-R", "UWAVE-SD",
                                           "MVLUE-R", "MVLUE-SD", "RMSDF"),
                             conf.nsigma  =  3, limits = NULL, plot = FALSE, ...)
#.........................................................................
  {
  
  if (!is.numeric(std.dev))
    std.dev <- match.arg(std.dev)

  obj<-qcd(data= x, var.index = var.index, sample.index = sample.index,
       covar.index = covar.index, covar.names = covar.names,
       data.name = data.name, sizes = sizes)

  result<-qcs.xbar.qcd(x = obj,  center = center, std.dev = std.dev,
                       conf.nsigma = conf.nsigma, 
                       limits = limits, plot = plot)

  return(result)
} # qcs.xbar.default
#.........................................................................

##' @rdname  qcs.xbar
##' @method qcs.xbar qcd
##' @inheritParams qcs.xbar.default
##' @export
##' 
qcs.xbar.qcd <- function(x, center = NULL,
                         std.dev  =  c("UWAVE-R", "UWAVE-SD",
                                       "MVLUE-R", "MVLUE-SD", "RMSDF"),
                         conf.nsigma  =  3, limits = NULL, plot = FALSE, ...) 
#.........................................................................  
{

  if (!is.numeric(std.dev))
    std.dev <- match.arg(std.dev)
  
  if(is.null(x) || !inherits(x, "qcd"))
    stop("data must be an objects of class (or extending) 'qcd'")
  
  sizes <- x$sizes
  type.data <- "continuous"
  
  qcs<-qcs(x$x, x$sample, sizes, type  =  "xbar",
           center, std.dev, conf.nsigma, limits, type.data)
  
  center <- qcs$center
  xbar <- qcs$statistics
  std.dev <- qcs$std.dev
  sizes <- qcs$sizes
  limits <- qcs$limits
  violations <- qcs$violations
  
  statistics <- data.frame(xbar)
  m <- length(x)
  sample <- x$sample

  if (m > 3) {
    new.x <- x[, -c(1, 2, length(x))]
    cov <- apply(new.x, 2, function(x) unlist(lapply(split(x, sample), unique)))
    statistics <- data.frame(xbar, cov)
  }
  
  row.names(statistics) <- unique(x$sample)
  data.name <- attr(x, "data.name")
  result <- list(qcd  =  x, type  =  "xbar", statistics  =  statistics,
                 center  =  center, std.dev  =  std.dev,
                 limits  =  limits, conf.nsigma  =  conf.nsigma,
                 sizes  =  sizes, data.name  =  data.name,
                 violations  =  violations)
  
  oldClass(result) <- c("qcs.xbar", "qcs")
  
  if(plot) plot(result, ...)
  
  return(result)
  #.........................................................................
} # qcs.xbar.qcd
#.........................................................................

