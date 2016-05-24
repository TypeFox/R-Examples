# private function
# Estimates an appropiate number of boxes for using in the box assisted algorithm
# using the time series being analyzed and the radius of the grid.
# @details
# The estimation takes the difference between the maximum and the minimum of the time
# series and divides it by the radius of the grid. If the number of boxes is too large
# (over a kMaximumNumberBoxes), the number of boxes is fixed to kMaximumNumberBoxes.
# If the number of boxes is too small (below kMinimumNumberBoxes), the number of boxes
# is fixed to kMinimumNumberBoxes.
# @param time.series the original time.series from which the surrogate data is generated
# @param radius width of each of the  rectangles that will be used to build the grid
# in the box assisted algorithm
# @return the number of boxes to used in the box assisted algorithm
# @author Constantino A. Garcia
estimateNumberBoxes = function(time.series, radius){
   kMinimumNumberBoxes = 10
   kMaximumNumberBoxes = 500
   number.boxes = as.integer( ( max(time.series)-min(time.series) )/radius  )
   if (number.boxes > kMaximumNumberBoxes) number.boxes = kMaximumNumberBoxes
   if (number.boxes < kMinimumNumberBoxes) number.boxes = kMinimumNumberBoxes
   
   return (number.boxes)
}

# private function
# checks if takens' vectors v1 and v2 are neighbours
# the neighbourhood is determined using the max norm and an radius radious.
# embedding.dim is not strictly necessary but it is used to accelerate the computations
isNeighbour=function(v1,v2,embedding.dim,radius){
  for (i in 1:embedding.dim){
    if ( abs(v1[[i]]-v2[[i]]) >=radius) return (FALSE);
  }
  return (TRUE);
}


################################################################################
#' Build the Takens' vectors 
#' @description
#' This function builds the Takens' vectors from a given time series. The set
#' of Takens' vectors is the result of embedding the time series in a m-dimensional
#' space. That is, the \eqn{n^{th}} Takens' vector is defined as 
#' \deqn{T[n]=\{time.series[n], time.series[n+ timeLag],...time.series[n+m*timeLag]\}.}
#' Taken's theorem states that we can then reconstruct an equivalent dynamical 
#' system to the original one (the 
#' dynamical system that generated the observed time series) by using the Takens' vectors.

#' @param time.series The original time series.
#' @param embedding.dim Integer denoting the dimension in which we shall embed the time.series.
#' @param time.lag Integer denoting the number of time steps that will be use to construct the 
#' Takens' vectors.
#' @return A matrix containing the Takens' vectors (one per row). The resulting
#' matrix also contains information about the time lag and the embedding
#' dimension used (as attributes).
#' @references H. Kantz  and T. Schreiber: Nonlinear Time series Analysis (Cambridge university press)
#' @author Constantino A. Garcia and Gunther Sawitzki.
#' @examples 
#' \dontrun{
#'# Build the Takens vector for the Henon map using the x-coordinate time series
#' h = henon(n.sample=  3000,n.transient= 100, a = 1.4, b = 0.3, 
#' start = c(0.73954883, 0.04772637), do.plot = FALSE)
#' takens = buildTakens(h$x,embedding.dim=2,time.lag=1)
#' # using the x-coordinate time series we are able to reconstruct
#' # the state space of the Henon map
#' plot(takens)}
#' @export buildTakens
buildTakens=function (time.series, embedding.dim, time.lag) 
{
  N = length(time.series)
  #vector of increments
  maxjump = (embedding.dim - 1) * time.lag
  jumpsvect = seq(0, maxjump, time.lag)
  #matrix that will store the takens' vectos. One vector per row
  numelem = N - maxjump
  takens = matrix(nrow = numelem, ncol = embedding.dim)
  #build takens' vectors
  for (i in 1:numelem) {
    takens[i, 1:embedding.dim] = time.series[jumpsvect + i]
  }
  # 
  #class(takens) = "takens"
  id=deparse(substitute(time.series))
  attr(takens,"embedding.dim") = embedding.dim
  attr(takens,"time.lag") = time.lag
  attr(takens,"id") = id
  return(takens)
}


#' Estimate an appropiate time lag for the Takens' vectors
#' @description
#' Given a time series (time.series), an embedding dimension (m) and a 
#' time lag (timeLag), the \eqn{n^{th}} 
#' Takens' vector is defined as 
#' \deqn{T[n]=\{time.series[n], time.series[n+ timeLag],...time.series[n+m*timeLag]\}.}
#' This function estimates an appropiate time lag by using the autocorrelation function
#' or the average mutual information .
#' @details 
#' A basic criteria for estimating a proper time lag is based on the following reasoning:
#' if the time lag used to build the Takens' vectors is too small, the coordinates will
#' be too highly temporally correlated and the embedding will tend to cluster around 
#' the diagonal in the phase space. If the time lag is chosen too large, the resulting 
#' coordinates may be almost uncorrelated and the resulting embedding will be very complicated. 
#' Thus, the autocorrelation function can be used for  estimating an appropiate time lag of
#' a time series. However, it must be noted that the autocorrelation is a linear statistic,
#' and thus it does not take into account nonlinear dynamical correlations. To take into
#' account nonlinear correlations the average mutual information (AMI) can be used. 
#' Independently of the technique used to compute the correlation, the time lag can
#'  be selected in a variety of ways:   
#' \itemize{
#'    \item Select the time lag where the autocorrelation/AMI function decays to 0 
#'    (\emph{first.zero} selection method). This
#'    method is not appropriate for the AMI function, since it only takes positive values.
#'    \item Select the time lag where the autocorrelation/AMI function decays to 
#'    1/e of its value at zero (\emph{first.e.decay} selection method).
#'    \item Select the time lag where the autocorrelation/AMI function reaches 
#'    its first minimum (\emph{first.minimum} selection method).
#'    \item Select the time lag where the autocorrelation/AMI function decays to
#'     the value specified by the user (\emph{first.value} selection method and 
#'     \emph{value} parameter).
#' }   
#' @param time.series The original time series.
#' @param technique The technique that we shall use to estimate the time lag (see the 
#' Details section). Allowed values are \emph{"acf"} and \emph{"ami"}.
#' @param selection.method Method used for selecting a concrete time lag. 
#' Available methods are \emph{"first.zero"}, \emph{"first.e.decay"} (default),
#'  \emph{"first.minimum"} and \emph{"first.value"}. 
#' @param value Numeric value indicating the value that the autocorrelation/AMI 
#' function must cross in order to select the time lag. It is used only with
#'  the "first.value" selection method.
#' @param lag.max Maximum lag at which to calculate the acf/AMI. 
#' @param do.plot Logical value. If TRUE (default value), a plot of the 
#' autocorrelation/AMI function is shown.
#' @param main A title for the plot.
#' @param ... Additional parameters for the \emph{acf} or the 
#' \emph{mutualInformation}.
#' @return The estimated time lag.
#' @examples
#' \dontrun{
#' sx = sinaiMap(a=0.3,n.sample=5000,start=c(0.23489,0.8923),do.plot=FALSE)$x
#' timeLag(sx, technique="ami",  
#'         n.partitions = 20, units = "Bits")
#' timeLag(sx, technique="acf") }
#' @note If the autocorrelation/AMI function does not cross the specifiged value, 
#' an error is thrown. This may be solved by increasing the \emph{lag.max} or 
#' selecting a higher value to which the autocorrelation/AMI function may decay.
#' @seealso \code{\link{mutualInformation}}
#' @references H. Kantz  and T. Schreiber: Nonlinear Time series Analysis (Cambridge university press)
#' @author Constantino A. Garcia
#' @export timeLag
timeLag = function (time.series, technique=c("acf","ami"), 
                    selection.method=c("first.e.decay", "first.zero", 
                                       "first.minimum", "first.value"),
                    value = 1/exp(1), lag.max=NULL, do.plot=TRUE,
                    main=NULL,...){
  
  technique = match.arg(technique)
  if (is.null(main)){
    main = switch(technique,
                  acf = "Autocorrelation function",
                  ami = "Average Mutual Information (AMI)")
  }
  if( is.null(lag.max) && (technique == "acf") ){
    
    lag.max = max(20, length(time.series)/2)
    # if technique=="ami", the mutualInformation function will handle the null value
  }
  fx = switch(technique, 
              ami = mutualInformation(time.series,lag.max = lag.max,
                                      do.plot=do.plot,main=main,... )$mutual.information,
              acf = stats::acf(time.series,lag.max=lag.max,plot=do.plot,
                               type="correlation",main=main,...)$acf
              )
  
  
  selection.method = match.arg(selection.method)
  lag = switch(match.arg(selection.method),
               first.e.decay = get.position.decays(fx,fx[1]/exp(1),technique),
               first.zero = get.position.decays(fx,0,technique),
               first.minimum = get.first.minimum(fx,technique),
               first.value = get.position.decays(fx,value,technique) )
  if (do.plot){
    if (selection.method == "first.e.decay"){
      abline(h=fx[1]/exp(1),lty=3,col=3,lwd=2)
    }else if(selection.method == "first.value"){
      abline(h=value,lty=3,col=3,lwd=2)
    }
  } 
  
  lag
}


# internal function. It will calculate the first time lag where
# the vector x decays to a given value
get.position.decays=function(x, value, method ){
  cross.position <- which(x <= value)
  # if x does not cross the value specified
  # by the user, a warning is given
  if (length(cross.position) == 0){
    stop(paste("The ",method," function does not cross ", value,
               ". Choose another \"cross\" value!\n",sep=""))
  }else{   
    cross.position = cross.position[[1]]
  }
  # If possible: convert from positions to time lags by substracting 1
  if ((cross.position) > 1){
    cross.position = cross.position-1
  }
  return (cross.position)
}

# internal function. It will calculate the first time lag where
# the x vector has a minimum
get.first.minimum=function(x,method){
  derivative = diff(x)
  derivative = ifelse( derivative < 0, -1, derivative)
  derivative = ifelse( derivative == 0, 0, derivative)
  derivative = ifelse( derivative > 0, 1, derivative)
  minimums = which(diff(derivative) >= 1) 
  # if the autocorrelation function is monotonically decreasing, an error is given
  if (length(minimums) == 0) {
    stop(paste("The ",method," function does not have a minimum\n",sep=""))
  }else{
    #take into account that we have loss some samples after diff
    
    first.minimum = minimums[[1]] + 1
  } 
  # convert from positions to time lags by substracting 1
  return (first.minimum-1)
}
