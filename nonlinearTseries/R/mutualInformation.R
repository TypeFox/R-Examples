#' Average Mutual Information (AMI)
#' @description
#' Functions for estimating the Average Mutual Information (AMI) of a time
#' series.
#' @details 
#' The Average Mutual Information (AMI) measures how much one random variable tells us about 
#' another. In the context of time series analysis, AMI
#' helps to quantify the amount of knowledge gained about the value
#' of \eqn{x(t+\tau)}{x(t+tau)} when observing \eqn{x(t)}.
#' 
#' To measure the AMI iof a time series, we create a histogram of the data
#' using bins. Let \eqn{p_i}{Pi} the probability that the signal has a 
#' value inside the ith bin, and let  \eqn{p_{ij}(\tau)}{Pij(tau)} be
#' the probability that \eqn{x(t)} is in bin i ans \eqn{x(t+\tau)}{x(t+tau)}
#' is in bin j. Then, AMI for time delay \eqn{\tau}{tau} is defined as
#' 
#' \deqn{AMI(\tau) = \sum_{i,j} p_{ij} log(\frac{p_{ij}}{p_i p_j})}{
#' AMI(tau) = sum( Pij log( Pij / (Pi*Pj) ) ) }
#' 
#' Depending on the base of the logarithm used to define AMI, the AMI
#' is measured in bits (base 2, also called shannons), nats (base e) or
#'  bans (base 10, also called hartleys).
#' @param time.series The observed time series.
#' @param lag.max Largest lag at which to calculate the AMI.
#' @param n.partitions Number of bins used to compute the probability distribution
#' of the time series.
#' @param units The units for the mutual information. Allowed units are
#' "Nats", "Bits" or "Bans" (somethings called Hartleys). Default is "Nats".
#' @param do.plot Logical value. If TRUE, the AMI is plotted
#' @param ... Further arguments for the plotting function.
#' @return  A \emph{mutualInf} object that consist of a list containing all 
#' the relevant information of the AMI computation: 
#' \emph{time.lag}, \emph{mutual.information}, \emph{units} and \emph{n.partitions}. 
#' @references H. Kantz  and T. Schreiber: Nonlinear Time series Analysis (Cambridge university press)
#' H. Abarbanel: Analysis of observed chaotic data (Springer, 1996).
#' @author Constantino A. Garcia
#' @examples 
#' \dontrun{
#' sx = sinaiMap(a=0.3,n.sample=5000,start=c(0.23489,0.8923),do.plot=FALSE)$x
#' mutinf = mutualInformation(sx, n.partitions = 20, units = "Bits") }
#' @seealso \code{\link{timeLag}}
#' @rdname mutualInformation 
#' @export mutualInformation
#' @exportClass mutualInf
#' @useDynLib nonlinearTseries
#' @import Rcpp
mutualInformation = function(time.series,lag.max = NULL, 
                              n.partitions = NULL, 
                              units = c("Nats","Bits","Bans"),
                              do.plot=TRUE,
                              ...){
  units = match.arg(units)
  base = switch(units, 
                Nats = exp(1),
                Bits = 2,
                Bans = 10)
  
  if(is.null(lag.max)){
    lag.max = max(20, sqrt(length(time.series)) )
  }  
  if (is.null(n.partitions)){
    n.partitions = max(floor(length(time.series)^(1/3)),2)
  }
  #call C
  mutinf = .Call("nonlinearTseries_mutualInformation", as.numeric(time.series),
                  as.integer(lag.max),
                  as.integer(n.partitions),base=as.numeric(base),
                  PACKAGE="nonlinearTseries" )
  mutinf = mutualInf(0:lag.max,mutinf / log(base), units, n.partitions)
  
  if(do.plot){
    plot(mutinf)
  }
  
  mutinf
}

#' @param x A \emph{mutualInf} object.
#' @param main Title for the plot.
#' @param xlab Title for the x axis.
#' @param ylab Title for the y axis.
#' @param type Type of plot to be drawn.
#' @rdname mutualInformation
#' @method plot mutualInf
#' @export
plot.mutualInf = function(x,main="Average Mutual Information (AMI)",
                          xlab="Time lag", ylab=NULL, type="h",
                          ...){
  if (is.null(ylab)){
    ylab = paste("AMI in", x$units)
  }
  x.len = length(x$mutual.information)
  plot(x$time.lag, x$mutual.information,type=type, main=main,
       xlab=xlab,ylab=ylab,...)
}


#' @rdname mutualInformation
#' @method as.numeric mutualInf
#' @export
as.numeric.mutualInf = function(x,...){
  x$mutual.information
}



# Extract parts of the AMI object
#' @param i Indices specifying elements to extract.
#' @rdname mutualInformation
#' @export
`[.mutualInf` = function(x,i){
    mutualInf(x$time.lag[i],x$mutual.information[i], x$units, x$n.partitions)
}


#  Extract parts of the AMI object
#' @rdname mutualInformation
#' @export
`[[.mutualInf` = function(x,i){
  mutualInf(x$time.lag[i], x$mutual.information[i], x$units, x$n.partitions)
}


mutualInf = function(tim, mutinf, units, n.partitions){
  mutinf = list(time.lag = tim, mutual.information = mutinf, units=units, 
                n.partitions = n.partitions)
  class(mutinf) = "mutualInf"
  mutinf
}

tsHist <- function(x,tlag=1,npartitions=16){
  .Call("nonlinearTseries_tsHistogram",as.numeric(x),as.integer(tlag),
        as.integer(npartitions),PACKAGE="nonlinearTseries")
}

