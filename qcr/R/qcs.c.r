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
# c chart
#-------------------------------------------------------------------------
##' Function to plot Shewhart c chart
##'
##' This function is used to compute statistics required by the c chart.
##'
##' @param x   an R object (used to select the method). See details.
##' @param ... arguments passed to or from methods.
##' @export
##' @examples
##' library(qcr)
##' data(circuit)
##' attach(circuit)
##' str(circuit)
##' datos <- circuit
##' datos$sample <- 1:length(datos$x)
##' str(datos)
##' sizes <- datos[,2]
##'
##' datos.qcd <- qcd(data = datos, var.index = 1,sample.index = 2,
##'                 sizes = size, type.data = "atributte")
##' res.qcs <- qcs.c(datos.qcd)
##' summary(res.qcs)
##' plot(res.qcs)
##'


qcs.c <- function(x, ...) {
  UseMethod("qcs.c")
}

##' @rdname qcs.c
##' @method qcs.c default
##' @inheritParams qcd
##' @param center a value specifying the center of group statistics or the
##' ''target'' value of the process.
##' @param conf.nsigma  a numeric value used to compute control limits, specifying the
##' number of standard deviations (if \code{conf.nsigma} > 1) or the confidence level (if 0
##' < \code{conf.nsigma} < 1).
##' @param limits a two-values vector specifying control limits.
##' @param plot a logical value indicating should be plotted.
##' @export
qcs.c.default <- function(x, var.index  =  1, sample.index  =  2,
                             covar.index  =  NULL, covar.names  =  NULL,
                             data.name = NULL,
                             sizes  =  NULL,
                             center = NULL,
                          conf.nsigma  =  3, limits = NULL, plot = FALSE, ...)
  {
  if (is.null(sizes)) 
    stop("sample sizes must be given for a attribute variable")


  obj<-qcd(data = x, var.index = var.index, sample.index = sample.index,
       covar.index = covar.index, covar.names = covar.names,
       data.name = data.name, sizes = sizes, type.data = "atributte")

  result<-qcs.c.qcd(x = obj, center = center, conf.nsigma = conf.nsigma,
                    limits = limits, plot = plot)

  return(result)
}


##' @rdname  qcs.c
##' @method qcs.c qcd
##' @inheritParams qcs.c.default
##' @export
qcs.c.qcd <- function(x, center = NULL,
                         conf.nsigma  =  3, limits = NULL, plot = FALSE, ...) {
  #.........................................................................
  if(is.null(x) || !inherits(x, "qcd"))
    stop("data must be an objects of class (or extending) 'qcd'")
  sizes <- x$sizes
  type.data <- "atributte"
  
  qcs<-qcs(x = x$x, sample.index = x$sample, sizes = sizes, type  =  "c",
            center = center, 
           conf.nsigma = conf.nsigma, limits = limits, type.data = type.data)
  
  center <- qcs$center
  c <- qcs$statistics
  std.dev <- qcs$std.dev
  sizes <- qcs$sizes
  limits <- qcs$limits
  violations <- qcs$violations
  
  statistics <- data.frame(c)
  m <- length(x)
  sample <- x$sample

  if (m > 3) {
    new.x <- x[, -c(1, 2, length(x))]
    cov <- apply(new.x, 2, function(x) unlist(lapply(split(x, sample), unique)))
    statistics <- data.frame(c, cov)
  }
  
  row.names(statistics) <- unique(x$sample)
  data.name <- attr(x, "data.name")
  result <- list(qcd  =  x, type  =  "c", statistics  =  statistics,
                 center  =  center, std.dev  =  std.dev,
                 limits  =  limits, conf.nsigma  =  conf.nsigma,
                 sizes  =  sizes, data.name  =  data.name,
                 violations  =  violations)
  
  oldClass(result) <- c("qcs.c", "qcs")
  
  if(plot) plot(result, ...)
  
  return(result)
  #.........................................................................
} # qcs.u.qcd