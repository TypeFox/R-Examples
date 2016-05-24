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
  # S chart
  #-------------------------------------------------------------------------
  ##' Function to plot Shewhart S chart
  ##'
  ##' This function is used to compute statistics required by the S chart.
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
  ##'data(pistonrings)
  ##'str(pistonrings)
  ##'pistonrings.qcd<-qcd(pistonrings)
  ##'
  ##'class(pistonrings.qcd)
  ##'
  ##'res.qcs <- qcs.S(pistonrings.qcd)
  ##'class(res.qcs)
  ##'plot(res.qcs,title="Control Chart S for pistonrings")
  ##'summary(res.qcs)  
  ##'
  qcs.S <- function(x, ...) {
    UseMethod("qcs.S")
  }

  ##' @rdname qcs.S
  ##' @method qcs.S default
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
  ##' @details
  ##' In the default method \code{qcs.S.default} parameter \code{x} is a matrix
  ##' or data-frame where it should contain data, index sample and, optionally, covariate(s).
  ##' @seealso
  ##' \code{\link{qcs}}, \code{\link{qcd}}
  ##' @export
  qcs.S.default <- function(x, var.index  =  1, sample.index  =  2,
                            covar.index  =  NULL, covar.names  =  NULL,
                            data.name = NULL,
                            sizes  =  NULL,
                            center = NULL,
                            std.dev  =  c("UWAVE-SD", "MVLUE-SD", "RMSDF"),
                            conf.nsigma  =  3, limits = NULL, plot = FALSE, ...)
  {
    if (!is.numeric(std.dev))
      std.dev <- match.arg(std.dev)

    obj<-qcd(data= x, var.index = var.index, sample.index = sample.index,
             covar.index = covar.index, covar.names = covar.names,
             data.name = data.name, sizes = sizes)
    
    result<-qcs.S.qcd(x = obj,  center = center, std.dev = std.dev,
                         conf.nsigma = conf.nsigma, 
                         limits = limits, plot = plot)   
    return(result)
  }



  ##' @rdname  qcs.S
  ##' @method qcs.S qcd
  # @inheritParams qcs.S.default
  ##' @export
  qcs.S.qcd <- function(x, center = NULL,
                        std.dev = c("UWAVE-SD", "MVLUE-SD", "RMSDF"),
                        conf.nsigma  =  3, limits = NULL, plot = FALSE, ...) {
    #.........................................................................
    if (!is.numeric(std.dev))
      std.dev <- match.arg(std.dev)


    if(is.null(x) || !inherits(x, "qcd"))
      stop("data must be an objects of class (or extending) 'qcd'")

    sizes <- x$sizes
    type.data <- "continuous"
    

    qcs<-qcs(x$x, x$sample, sizes, type  =  "S",
             center, std.dev, conf.nsigma, limits, type.data)

    center <- qcs$center
    S <- qcs$statistics
    std.dev <- qcs$std.dev
    sizes <- qcs$sizes
    limits <- qcs$limits
    violations <- qcs$violations

    statistics <- data.frame(S)
    m <- length(x)
    sample <- x$sample
    if (m > 3) {
      new.x <- x[, -c(1, 2, length(x))]
      cov <- apply(new.x, 2, function(x) unlist(lapply(split(x, sample), unique)))
      statistics <- data.frame(S, cov)
    }

    row.names(statistics) <- unique(x$sample)
    data.name <- attr(x, "data.name")
    result <- list(qcd  =  x, type  =  "S", statistics  =  statistics,
                   center  =  center, std.dev  =  std.dev,
                   limits  =  limits, conf.nsigma  =  conf.nsigma,
                   sizes  =  sizes, data.name  =  data.name,
                   violations  =  violations)

    oldClass(result) <- c("qcs.S", "qcs")

    if(plot) plot(result, ...)

    return(result)
    #.........................................................................
  } # qcs.S.qcd
