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
  # one chart
  #-------------------------------------------------------------------------
  ##' Function to plot Shewhart xbar.one chart
  ##' 
  ##' This function is used to compute statistics required by the xbar.one chart.
  ##'   
  ##' @param x   an R object (used to select the method). See details.
  ##' @param ... arguments passed to or from methods.
  ##' @export
  ##' @examples
  ##' 
  ##' ##
  ##' ##  Continuous data 
  ##' ##
  ##'library(qcr)
  ##' x <- c(33.75, 33.05, 34, 33.81, 33.46, 34.02, 33.68, 33.27, 33.49, 33.20,
  ##'       33.62, 33.00, 33.54, 33.12, 33.84)
  ##' 
  ##' sample <- 1:length(x)
  ##' datos <- data.frame(x,sample)
  ##' datos.qcd <- qcd(datos)
  ##'
  ##' res.qcs <- qcs.one(datos.qcd)
  ##' class(res.qcs)
  ##' summary(res.qcs)  
  ##'  plot(res.qcs, title = "Control Chart Xbar.one for pistonrings")  
  ##'
  
  qcs.one <- function(x, ...) {
    UseMethod("qcs.one")
  }
  
  ##' @rdname  qcs.one
  ##' @method qcs.one default
  ##' @inheritParams qcd
  ##' @param sizes  optional. A value or a vector of values specifying the sample sizes
  ##' associated with each group. For continuous data the sample sizes are obtained counting the non-\code{NA} elements of
  ##' the sample.index vector. For attribute
  ##' variable the argument sizes is required.
  ##' @param center a value specifying the center of group statistics or the
  ##' ''target'' value of the process.
  ##' @param std.dev  a value or an available method specifying the within-group standard
  ##' deviation(s) of the process. Several methods are available for estimating the
  ##' standard deviation in case of a continuous process variable.
  ##' @param k number of successive pairs of observations for computing the
  ##' standard deviation based on moving ranges of k points.
  ##' @param conf.nsigma  a numeric value used to compute control limits, specifying the
  ##' number of standard deviations (if \code{conf.nsigma} > 1) or the confidence level (if 0
  ##' < \code{conf.nsigma} < 1).
  ##' @param limits a two-values vector specifying control limits.
  ##' @param plot a logical value indicating should be plotted.
  ##' @details
  ##' In the default method \code{qcs.one.default} parameter \code{x} is a matrix 
  ##' or data-frame where it should contain data, index sample and, optionally, covariate(s).  
  ##' @export
qcs.one.default <- function(x, var.index  =  1, sample.index  =  2, 
                                   covar.index  =  NULL, covar.names  =  NULL,
                                   data.name = NULL, 
                                   sizes = NULL,
                                   center = NULL, 
                                   std.dev = c("MR", "SD"), k = 2,
                                   conf.nsigma  =  3, 
                                   limits = NULL, plot = FALSE, ...) 
  { 
    std.dev <- match.arg(std.dev)
    
    obj <- qcd(data = x, var.index = var.index, sample.index = sample.index, 
             covar.index = covar.index, covar.names = covar.names,
             data.name = data.name, sizes = sizes)  
    
    result <- qcs.one.qcd(x = obj, center = center,  std.dev = std.dev, 
                               k = k , conf.nsigma = conf.nsigma, 
                               limits = limits, plot = plot)
    
    return(result)
  }
  


  ##' @rdname  qcs.one
  ##' @method qcs.one qcd
  ##' @inheritParams qcs.one.default
  ##' @export
  qcs.one.qcd <- function(x, center = NULL, 
                               std.dev = c("MR", "SD"), k = 2,
                               conf.nsigma  =  3, 
                               limits = NULL, plot = FALSE, ...) {
    #.........................................................................
    if (!is.numeric(std.dev))
      std.dev <- match.arg(std.dev)
    
  
    if(is.null(x) || !inherits(x, "qcd"))
      stop("data must be an objects of class (or extending) 'qcd'")
    type.data <- "continuous"
    std.dev <- sd.xbar.one (data = x$x, std.dev = std.dev, k = k)
    

    qcs <- qcs(x = x$x, sample.index = x$sample, type  =  "one", std.dev = std.dev,
             center = center, conf.nsigma = conf.nsigma, 
               limits = limits, 
               type.data = type.data)

    
    center <- qcs$center
    one <- qcs$statistics
    std.dev <- qcs$std.dev
    sizes <- qcs$sizes
    limits <- qcs$limits
    violations <- qcs$violations
    
    statistics <- data.frame(one)
    m <- length(x)
    sample <- x$sample 
    if (m > 3) {
      new.x <- x[, -c(1, 2, length(x))]
      cov <- apply(new.x, 2, function(x) unlist(lapply(split(x, sample), unique)))
      statistics <- data.frame(one, cov)
    }
    
    row.names(statistics) <- unique(x$sample)
    data.name <- attr(x, "data.name")
    result <- list(qcd  =  x, type  =  "one", statistics  =  statistics, 
                   center  =  center, std.dev  =  std.dev, 
                   limits  =  limits, conf.nsigma  =  conf.nsigma, 
                   sizes  =  sizes, data.name  =  data.name, 
                   violations  =  violations)
    
    oldClass(result) <- c("qcs.one", "qcs")
    
    if(plot) plot(result, ...) 
    
    return(result)
    #.........................................................................  
  } # qcs.one.qcd