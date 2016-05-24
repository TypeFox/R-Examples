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
#
#  Main function to create a 'qcs' object
#
##' Quality Control Statistics
##' 
##' Create an object of class 'qcs' to perform statistical quality control.
##' This object may then be used to plot Shewhart charts, Multivariate Control Charts,
##' and more.
##' 
##' 
##' @aliases qcs summary.qcs print.qcs
##' 
##' @param x  a vector containing observed data
##' @param sample.index a scalar with the column number corresponding the index each
##' group (sample).
##' @param sizes a value or a vector of values specifying the sample sizes
##' associated with each group. For continuous data the sample sizes are obtained counting the non-\code{NA} elements of
##' the sample.index vector. For \code{"p"}, \code{"np"} and \code{"u"} charts the argument
##' \code{sizes} is required.
##' @param type a character string specifying the group statistics to compute:
##' 
##' \tabular{lll}{ \tab Statistic charted \tab Chart description \cr
##' \code{"xbar"} \tab mean \tab means of a continuous process variable \cr
##' \code{"R"} \tab range \tab ranges of a continuous process variable \cr
##' \code{"S"} \tab standard deviation \tab standard deviations of a continuous
##' variable \cr \code{"one"} \tab mean \tab one-at-time data of a
##' continuous process variable \cr \code{"p"} \tab proportion \tab proportion
##' of nonconforming units \cr \code{"np"} \tab count \tab number of
##' nonconforming units \cr \code{"c"} \tab count \tab nonconformities per unit
##' \cr \code{"u"} \tab count \tab average nonconformities per unit \cr
##' \code{"g"} \tab count \tab number of non-events between events }
##' @param center a value specifying the center of group statistics or the
##' ''target'' value of the process.
##' @param std.dev  a value or an available method specifying the within-group standard
##' deviation(s) of the process. Several methods are available for estimating the
##' standard deviation in case of a continuous process variable.
##' @param conf.nsigma  a numeric value used to compute control limits, specifying the
##' number of standard deviations (if \code{conf.nsigma} > 1) or the confidence level (if 0
##' < \code{conf.nsigma} < 1).
##' @param limits a two-values vector specifying control limits.
##' @param type.data  a string specifying el type de data.
##' @param lambda the smoothing parameter \eqn{0 \le \lambda \le 1}{0 <= lambda
##' <= 1}
##' @param nsigma  a numeric value used to compute control limits, specifying the
##' number of standard deviations.
##' @param decision.interval A numeric value specifying the number of standard
##' errors of the summary statistics at which the cumulative sum is out of
##' control.
##' @param se.shift The amount of shift to detect in the process, measured in
##' standard errors of the summary statistics.
##' @return Returns an object of class 'qcs'.
##' @references Montgomery, D.C. (2000) \emph{Introduction to Statistical
##' Quality Control}, 4th ed. New York: John Wiley & Sons. \cr Wetherill, G.B.
##' and Brown, D.W. (1991) \emph{Statistical Process Control}. New York:
##' Chapman & Hall.
##' @export

#.........................................................................
qcs<-function(x, sample.index, sizes = NULL, type = c("xbar", "R", "S", "one", "p", "np", "c", "u", "ewma", "cusum"),
              center = NULL, std.dev, 
              conf.nsigma = 3, limits = NULL, type.data  = c("continuous","atributte","dependence"),
              lambda = 0.2, decision.interval = 5 ,
              se.shift = 1)              
#.........................................................................  
{

  type.data <- match.arg(type.data)
  type <- match.arg(type)

  
  x.qcc <- switch(type.data, 
                  "continuous" =  {qcs.continuous(x = x, sample.index = sample.index, sizes = sizes,
                                    type = type, center = center, 
                                    std.dev = std.dev, conf.nsigma = conf.nsigma,
                                    limits = limits)} ,
                  "atributte" =   {qcs.atributte(x = x, sample.index = NULL,
                                    sizes = sizes, type = type, center = center, 
                                    conf.nsigma = conf.nsigma, limits = limits)},
                  
                  "dependence" =  if (type == "ewma"){
                                     qcs.dependence(x = x, sample.index = sample.index, 
                                                    sizes =  sizes, type = type, center = center, 
                                                  std.dev = std.dev, nsigma = conf.nsigma, lambda = lambda)
                                   } else {
                                     qcs.dependence(x = x, sample.index = sample.index, sizes = sizes, 
                                                    type = type, center = center,
                                                    std.dev = std.dev, decision.interval = decision.interval,
                                                    se.shift = se.shift)})
  

    center <- x.qcc$center
    statistics <- x.qcc$statistics
    std.dev <- x.qcc$std.dev
    sizes <- x.qcc$sizes
    if (is.null(limits)) limits <- x.qcc$limits  
    violations <- x.qcc$violations

if (type.data == "dependence"){
    if (type == "ewma") {
          x <- x.qcc$x
          y <- x.qcc$y
          sigma <- x.qcc$sigma
          lambda <- x.qcc$lambda
          nsigma <- x.qcc$nsigma
        
          result <- list(statistics  =  statistics, center  =  center,
                         std.dev  =  std.dev, limits  =  limits, 
                         nsigma  =  nsigma, sizes  =  sizes,
                         violations  =  violations, x = x, y = y, lambda = lambda, 
                         sigma = sigma)
    } else {
      
      pos <- x.qcc$pos
      neg <- x.qcc$neg
      decision.interval <- x.qcc$decision.interval
      se.shift <- x.qcc$se.shift
      
      
      result <- list(statistics  =  statistics, center  =  center,
                     std.dev  =  std.dev, limits  =  limits, 
                     sizes  =  sizes,
                     violations  =  violations, pos = pos, 
                     neg = neg, decision.interval = decision.interval, se.shift = se.shift)      
    }
} else {
    result<-list(statistics  =  statistics, center  =  center,
               std.dev  =  std.dev, limits  =  limits, 
               sizes = sizes,
               violations  =  violations)
}

  oldClass(result) <- c("qcs")
  
  return(result)
} # qcs
#.........................................................................

##' @rdname  qcs

qcs.continuous<-function(x, sample.index, sizes = NULL, type = c("xbar", "R", "S", "one"), center = NULL, 
              std.dev, conf.nsigma = 3, limits = NULL)
  #.........................................................................  
{
  type = match.arg(type) 
  
  if (type != "one") {  
    x <- qcc.groups(x, sample.index)
    sizes <- table(sample.index)        
    if (is.null(center)) {
      if( is.null(limits)) {
        x.qcc<-qcc(data  =  x, sizes  =  sizes, type  =  type, 
                   std.dev  =  std.dev, confidence.level  =  conf.nsigma, plot  =  F )    
      } 
      else {
        x.qcc<-qcc(data  =  x, sizes  =  sizes, type  =  type, 
                   confidence.level  =  conf.nsigma,
                   plot  =  F, limits = limits)         
      }
    }          
    else {
      if (is.null(limits)){
        x.qcc<-qcc(data  =  x, sizes  =  sizes, center = center, type  =  type, 
                   std.dev  =  std.dev, confidence.level  =  conf.nsigma, plot  =  F)      
      } 
      else {
        x.qcc<-qcc(data  =  x, sizes  =  sizes, center = center, type  =  type, 
                   confidence.level  =  conf.nsigma, 
                   plot  =  F, limits = limits)
      }                 
    }     
  } 
  else {
    statistics <- as.vector(x)
    if (is.null(center)) 
      center <- mean(statistics)            
   
    if (is.null(limits)) {
      x.qcc <- qcc(data = statistics, type = "xbar.one", std.dev = std.dev, 
                   center = center, plot = F)
    } else {
      x.qcc <- qcc(data = statistics, type = "xbar.one", 
                   center = center,
                   limits = limits, plot = F)
    }
  }
  

  center <- x.qcc$center
  statistics <- x.qcc$statistics
  std.dev <- x.qcc$std.dev
  sizes <- x.qcc$sizes
  
  if (is.null(limits)) limits <- x.qcc$limits
  
  violations <- x.qcc$violations
  
  
  result<-list(statistics  =  statistics, center  =  center,  
               std.dev  =  std.dev, limits  =  limits, sizes = sizes,
               violations  =  violations)
  
  return( result)
  #   oldClass(result) <- c("qcs")
  #.........................................................................
} # qcs.continuous


##' @rdname  qcs

qcs.atributte<-function(x, sample.index = NULL, sizes = NULL, type = c("p", "np", "c", "u"), center = NULL, 
                         conf.nsigma = 3, limits = NULL)
  #.........................................................................  
{
  
  type = match.arg(type)
  
  if (is.null(sizes)) 
    stop("sample sizes must be given for a attribute variable")
  
  if (is.null(center)) {
    if( is.null(limits)) {
      x.qcc<-qcc(data  =  x, sizes  =  sizes, type  =  type, 
                 confidence.level  =  conf.nsigma, plot  =  F )    
    } 
    else {
      x.qcc<-qcc(data  =  x, sizes  =  sizes, type  =  type, 
                 confidence.level  =  conf.nsigma,
                 plot  =  F, limits = limits)         
    }
  } 
  else {
    if (is.null(limits)){
      x.qcc<-qcc(data  =  x, sizes  =  sizes, center = center, type  =  type, 
                 std.dev  =  std.dev, confidence.level  =  conf.nsigma, plot  =  F)      
    } 
    else {
      x.qcc<-qcc(data  =  x, sizes  =  sizes, center = center, type  =  type, 
                 confidence.level  =  conf.nsigma, 
                 plot  =  F, limits = limits)
    }                 
  }  
  

  
  center <- x.qcc$center
  statistics <- x.qcc$statistics
  std.dev <- x.qcc$std.dev
  sizes <- x.qcc$sizes
  
  if (is.null(limits)) limits <- x.qcc$limits
  
  violations <- x.qcc$violations
  
  
  result<-list(statistics  =  statistics, center  =  center,  
               std.dev  =  std.dev, limits  =  limits, sizes = sizes,
               violations  =  violations)
  
  return( result)

  #.........................................................................
} # qcs.atributte

##' @rdname  qcs

qcs.dependence<-function(x, sample.index = NULL, sizes = NULL, type = c("ewma","cusum"), center = NULL, 
                         std.dev, nsigma = 3, lambda = 0.2,  decision.interval = 5, se.shift = 1)
  #.........................................................................  
{
  
  type <- match.arg(type)
  if (type == "ewma") {    
    if (unique(sizes) == 1) {
        sizes <- 1
        if (is.null(center)) {
          x.qcc<-ewma(data  =  x, sizes  =  sizes, 
          std.dev  =  std.dev, nsigmas  =  nsigma, lambda = lambda , plot  =  F )    
        }
        else {
          x.qcc<-ewma(data  =  x, sizes  =  sizes, center = center, 
          std.dev  =  std.dev, nsigmas  =  nsigma, lambda = lambda , plot  =  F)         
          }
      } 
    else {
        x <- qcc.groups(x, sample.index)
        sizes <- table(sample.index)
        
        if (is.null(center)) {
          x.qcc<-ewma(data  =  x, sizes  =  sizes, 
          std.dev  =  std.dev, nsigmas  =  nsigma, lambda = lambda , plot  =  F )    
        }
        else {
          x.qcc<-ewma(data  =  x, sizes  =  sizes, center = center, 
          std.dev  =  std.dev, nsigmas  =  nsigma, lambda = lambda , plot  =  F)         
        }
        
      }
  }
  else {
    if (unique(sizes) == 1) {
      sizes <- 1
      if (is.null(center)) {
        x.qcc<-cusum(data  =  x, sizes  =  sizes, 
                    std.dev  =  std.dev, decision.interval  =  decision.interval, se.shift = se.shift , plot  =  F )    
      }
      else {
        x.qcc<-cusum(data  =  x, sizes  =  sizes, center = center, 
                    std.dev  =  std.dev, decision.interval  =  decision.interval, se.shift = se.shift , plot  =  F)         
      }
    } 
    else {
      x <- qcc.groups(x, sample.index)
      sizes <- table(sample.index)
      
      if (is.null(center)) {
        x.qcc<-cusum(data  =  x, sizes  =  sizes, 
                    std.dev  =  std.dev, decision.interval  =  decision.interval, se.shift = se.shift , plot  =  F )    
      }
      else {
        x.qcc<-cusum(data  =  x, sizes  =  sizes, center = center, 
                    std.dev  =  std.dev, decision.interval  =  decision.interval, se.shift = se.shift , plot  =  F)         
      }
      
    }
    
    
    
    
  }
  
  
  center <- x.qcc$center
  statistics <- x.qcc$statistics
  std.dev <- x.qcc$std.dev
  sizes <- x.qcc$sizes  
  violations <- x.qcc$violations

  if (type == "ewma") {
      limits <- x.qcc$limits
      xx <- x.qcc$x
      y  <- x.qcc$y
      sigma <- x.qcc$sigma
      lambda <- x.qcc$lambda
      nsigma  <-  x.qcc$nsigma
      result <- list(qcd  =  x, type  =  "ewma", statistics  =  statistics,
                     center  =  center, std.dev  =  std.dev,
                     limits  =  limits, nsigma  =  nsigma,
                     sizes  =  sizes,
                     violations  =  violations, x = xx, y = y, lambda = lambda, 
                     sigma = sigma)      
  } else  {  
      limits <- c(-x.qcc$decision.interval, x.qcc$decision.interval)
      pos <- x.qcc$pos
      neg <- x.qcc$neg
      decision.interval <- x.qcc$decision.interval
      se.shift <- x.qcc$se.shift
      result <- list(qcd  =  x, type  =  "cusum", statistics  =  statistics,
                     center  =  center, std.dev  =  std.dev,
                     limits  =  limits,
                     sizes  =  sizes,
                     violations  =  violations, 
                     pos = pos, neg = neg, decision.interval = decision.interval, 
                     se.shift = se.shift)
  }
  return( result)

  #.........................................................................
} # qcs.dependence
##' @export
print.qcs <- function(x, ...) str(x,1)
#.........................................................................
##' @export
summary.qcs <- function(object, ...)
#.........................................................................
{
  data.name <- object$data.name
  type <- object$type
  cat(paste(type, "chart for", data.name, "\n"))
  statistics <- object$statistics
  cat("\nSummary of group statistics:\n")
  print(summary(statistics))
  sizes <- object$sizes
  if(length(unique(sizes))==1)
    sizes <- sizes[1]
  if(length(sizes) == 1)
    cat("\nGroup sample size: ", format(sizes))
  else {
    cat("\nSummary of group sample sizes: ")
    tab <- table(object$sizes)
    print(matrix(c(as.numeric(names(tab)), tab), 
                 ncol = length(tab), byrow = TRUE, 
                 dimnames = list(c("  sizes", "  counts"),
                                 character(length(tab)))))
  }
  cat("\nNumber of groups: ", length(statistics[[1]]))
  center <- object$center
  cat("\nCenter of group statistics: ", center)
  sd <- object$std.dev
  cat("\nStandard deviation: ", sd, "\n")
  
  limits <- object$limits
  if (!is.null(limits)) 
  { cat("\nControl limits:", "\n") 
    print(limits)
  }
  
  beyond<-object$violations[[1]]
  violationg<-object$violations[[2]]
  
  if (length(object$violations[[1]])== 0){
    cat("\nNumber beyond limits: 0", "\n") 
  } 
  else {cat("\nBeyond limits of control:", "\n")
        print(object$statistics[beyond,])
  }
  
  if (length(object$violations[[2]])==0){
    cat("\nNumber violationg runs: 0", "\n") 
  } 
  else {cat("\nViolationg runs:", "\n")
        print(object$statistics[violationg,])
  }  
  invisible()
#.........................................................................
} # summary.qcs

# qcs.add function
#-------------------------------------------------------------------------
##' qcs.add Add a data.frame object with a qcs object
##' 
##' This function is used to join two objects of type data.frame and qcs.
##' 
##' @param x   Object type qcs
##' @export
##' 


qcs.add <- function(x, ...){
  UseMethod("qcs.add")
}

 
##' @rdname  qcs.add 
##' @method qcs.add default
##' @param value   Object type data.frame
##' @param var.index a scalar with the column number corresponding the observed data for
##' the variable (the variable quality).  Alternativelly can be a string with the
##' name of the quality variable.
##' @param sample.index a scalar with the column number corresponding the index each
##' group (sample).
##' @param covar.index  optional. A scalar or numeric vector with the column number(s)
##' corresponding to the covariate(s). Alternativelly can be a character vector with
##' the names of the covariates.
##' @param ...  arguments to be passed to or from methods.
##' @export 


qcs.add.default <- function(x, value, var.index = NULL,
                                 sample.index = NULL, covar.index = NULL, ...){

  
  if (!inherits(x, "qcs"))
    stop("object must be qcs")
  
  if (!is.matrix(value) & !is.data.frame(value))
    stop("object must be a matrix or data.frame")
  
  xx <- x$qcd
  center <- x$center
  std.dev <- x$std.dev
  limits <- x$limits
  type <- x$type
  
  if (length(xx)-1!=length(value))
    stop(" the objects must be the same length")
  
  if (is.null(var.index) & is.null(sample.index) & is.null(covar.index)) {
    yy <- value
  } else {
    yy <- value[c(var.index, sample.index, covar.index)]
  }
  
  sizes <- table(value[ ,2])
  yy$sizes <- sizes
  
  if (length(xx)==length(yy)){  
    names(yy) <- names(xx)    
    z <- rbind(xx,yy)
    n <- length(xx)
    
  }
  
  
  if (length(xx)>3){ 
    z.qcd<-qcd(data=z, covar.index = 3:length(z), 
               data.name = attr(xx,"data.name"), 
               type.data = attr(xx,"type.data"),
               sizes = z$sizes)
  } else {
    z.qcd<-qcd(data = z, data.name = attr(xx,"data.name"), 
               type.data = attr(xx,"type.data"),
               sizes = z$sizes)  
  }
  
  z.qcs <- switch(type, 
                  "xbar" = qcs.xbar.qcd(x = z.qcd, limits = limits),
                  "xbar.one" = qcs.one.qcd(x = z.qcd, limits = limits),
                  "R" = qcs.R.qcd(x = z.qcd, limits = limits),
                  "S" = qcs.S.qcd(x = z.qcd, limits = limits),
                  NULL)
  result <- z.qcs
}