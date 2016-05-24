########################### Determining embedding dimension ####################
#' Estimate the embedding dimension
#' @description
#' This function determines the minimum embedding dimension from a scalar time 
#' series using the algorithm proposed by L. Cao (see references).
#' @details
#' The Cao's algorithm uses 2 functions in order to estimate the embedding dimension
#' from a time series: the E1(d) and the E2(d) functions, where d denotes the dimension.
#' 
#' E1(d) stops changing when d is greater than or equal to the embedding dimension, staying close to 1.
#' On the other hand, E2(d) is used to distinguish deterministic signals from stochastic signals. For 
#' deterministic signals, there exist some d such that E2(d)!=1. For stochastic signals,
#' E2(d) is approximately 1 for all the values. 
#' 
#' This function uses the Arya and Mount's C++ ANN library for nearest neighbour search 
#' (For more information on the ANN library please visit \url{http://www.cs.umd.edu/~mount/ANN/}).
#' The R wrapper is a modified version of the RANN package code by Samuel E. Kemp and Gregory Jefferis.
#' @note
#' In the current version of the package, the automatic detection of stochastic 
#' signals has not been implemented yet.
#' @param time.series The original time series.
#' @param number.points Number of points from the time series that will be used to estimate
#' the embedding dimension. By default, all the points in the time series are used.
#' @param time.lag Time lag used to build the Takens' vectors needed to estimate the
#' embedding dimension (see \link{buildTakens}). Default: 1.
#' @param max.embedding.dim Maximum possible embedding dimension for the time series. Default: 15.
#' @param threshold Numerical value between 0 and 1. The embedding dimension is estimated
#' using the E1(d) function. E1(d) stops changing when d is greater than or equal to
#' embedding dimension, staying close to 1. This value establishes a threshold for 
#' considering that E1(d) is close to 1. Default: 0.95
#' @param max.relative.change Maximum relative change in E1(d) with respect to 
#' E1(d-1) in order to consider that the E1 function has been stabilized and it will
#' stop changing. Default: 0.01.
#' @param do.plot Logical value. If TRUE (default value), a plot of E1(d) and E2(d) is shown.
#' @param main Title for the plot.
#' @param xlab Title for the x axis.
#' @param ylab Title for the y axis.
#' @param xlim numeric vectors of length 2, giving the x coordinates range.
#' @param ylim numeric vectors of length 2, giving the y coordinates range.
#' @references 
#' Cao, L. Practical method for determining the minimum embedding dimension of a scalar time series. Physica D: Nonlinear Phenomena,
#' 110,1, pp. 43-50 (1997).
#' 
#' Arya S. and Mount D. M. (1993), Approximate nearest neighbor searching, Proc. 4th Ann. ACM-SIAM Symposium on Discrete Algorithms (SODA'93), 271-280.
#' 
#' Arya S., Mount D. M., Netanyahu N. S., Silverman R. and Wu A. Y (1998), An optimal algorithm for approximate nearest neighbor searching, Journal of the ACM, 45, 891-923.
#' @author Constantino A. Garcia
#' @examples 
#' \dontrun{
#' h = henon(do.plot=FALSE) 
#' dimension = estimateEmbeddingDim(h$x, time.lag=1, max.embedding.dim=6,
#'              threshold=0.9, do.plot=TRUE)
#'              }
#' @export estimateEmbeddingDim
estimateEmbeddingDim = function(time.series,  number.points = length(time.series), 
                                time.lag = 1,  max.embedding.dim = 15,  threshold = 0.95, 
                                max.relative.change = 0.10,
                                do.plot = TRUE,
                                main = "Computing the embedding dimension",
                                xlab="dimension (d)", ylab="E1(d) & E2(d)",
                                ylim = NULL, xlim = NULL){
  kSDFraction = 1e-6
  
  if (max.embedding.dim < 3) stop("max.embedding.dim should be greater that 2...\n")
  # normalize time series
  time.series = (time.series - mean(time.series,na.rm = T) ) / sd(time.series, na.rm = T)
  time.series.len = length(time.series)
  # add small quantity of noise to avoid finding identical phase space
  # points (something impossible in a pure chaotic signal). This "trick"
  # avoids failures when supplying periodic signals (i.e. a sine).
  # We use the IQR as an estimation of the time.series dispersion to
  # avoid the effect of outliers
  time.series = time.series +
    rnorm(time.series.len, 
          sd = IQR(time.series,na.rm = T) * kSDFraction) 
  data = time.series[(time.series.len/2-number.points/2+1):(time.series.len/2+number.points/2)]
  #if no d verifies E1(d) >= threshold,  then we shall return 0
  embedding.dim = 0 
  # First iteration: get E1(1) and E2(1)
  E.parameters = getCaoParameters(data,  1,  time.lag)
  E.parameters.next.dim = getCaoParameters(data,  2,  time.lag)
  E.vector = c(E.parameters$E,  E.parameters.next.dim$E)
  E.star.vector = c(E.parameters$E.star,  E.parameters.next.dim$E.star)
  E1.vector = c(E.vector[[2]]/E.vector[[1]])
  E2.vector = c(E.star.vector[[2]]/E.star.vector[[1]])
  # compute from d = 3 to d = max.embedding.dim
  for (dimension in 3:max.embedding.dim){
    #compute E parameters, E1 and E2
    E.parameters = getCaoParameters(data, dimension, time.lag)
    E.vector[[dimension]] = E.parameters$E
    E.star.vector[[dimension]] = E.parameters$E.star
    E1.vector[[dimension-1]] = E.vector[[dimension]]/E.vector[[dimension-1]]
    E2.vector[[dimension-1]] = E.star.vector[[dimension]]/E.star.vector[[dimension-1]]
    # Error for dimension - 2
    relative.error = abs(E1.vector[[dimension-1]]-E1.vector[[dimension-2]])/(E1.vector[[dimension-2]])
    #compute if E1(d)>=threshold...If it is the first time it happens(embedding.dim==0), store
    # the dimension
    if ((embedding.dim==0) && (!is.na(E1.vector[[dimension - 2]])) &&
          (E1.vector[[dimension-2]]>=threshold) &&
          (relative.error < max.relative.change )){
      embedding.dim = dimension - 2
    }
  }
  #plot graphics
  if (do.plot){
    if (is.null(ylim)){
      ylim = c(0,max(E1.vector,E2.vector))
    }
    if (is.null(xlim)){
      xlim = c(1,length(E1.vector))
    }
    plot(1:length(E1.vector), E1.vector, 'b',lty=1,pch=1,
         main=main,ylab=ylab,xlab=xlab, cex = 1, ylim = ylim, xlim=xlim)
    lines(1:length(E2.vector), E2.vector, 'b',lty=2,pch=2, col = 2, cex = 1)
    abline(h = c(1, threshold), col = 3, lwd=2,lty = c(3, 3))
    legend("bottomright", col = 1:3, lty = 1:3, lwd = c(2.5, 2.5, 2.5),
           pch=c(1,2,NA),bty="n",
           legend = c("E1(d)", "E2(d)", "limits for E1(d)"))
    
  }
  return (embedding.dim)
  
}


# private function
# auxiliar function to compute E,  E1 and E2 based on the 
# L.Cao article: Practical method for determining the minimum embedding dimension of a scalar time series.
getCaoParameters = function(data, m, time.lag){
  # theshold for considering that two vectors are at distance 0
  kZero = 1e-15
  #construct takens vectors of dimensions m and m+1
  takens = buildTakens(data,  m,  time.lag)
  takens.next.dimension = buildTakens(data,  m+1,  time.lag)
  #get closest neigh
  max.iter = nrow(takens.next.dimension)
  if (m==1){
    nearest.neigh = nn.search(data = matrix(data[1:max.iter],ncol=1), query = matrix(data[1:max.iter],ncol=1), k=2, eps=0)
  }else{
    nearest.neigh = nn.search(data = takens[1:max.iter,], query= takens[1:max.iter,], k=2, eps=0)
  } 
  # the a(i, d) parameter from the Cao's article (Equation 1) will be call here
  # min.dist.ratio. On the other hand,  the expression inside the summatory in 
  # equation 4 will be called stochastic.parameter
  min.dist.ratio = c()
  stochastic.parameter = c()
  #computing...
  for (takens.position in 1:max.iter){
    # get closest neighbour (avoid picking the same vector with the 2 index)
    closest.neigh = nearest.neigh$nn.idx[takens.position,2] 
    if (nearest.neigh$nn.dists[takens.position, 2] < kZero){
      # We found equal points in phase space... assing NA     
      min.dist.ratio[[takens.position]] = NA
    }else{
      numerator = as.numeric(dist(rbind(takens.next.dimension[takens.position, ], takens.next.dimension[closest.neigh, ]),  method = "maximum"))
      min.dist.ratio[[takens.position]] = numerator/nearest.neigh$nn.dists[takens.position, 2]
    }
    
    stochastic.parameter[[takens.position]] = abs(data[[takens.position+m*time.lag]]-data[[closest.neigh+m*time.lag]])
  }

  return (list(E = mean(min.dist.ratio,na.rm=TRUE), E.star = mean(stochastic.parameter,na.rm=TRUE)))
}
