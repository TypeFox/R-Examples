#' @include generics.R
#' @include Class-LagOperator.R

NULL

################################################################################
#' Class to calculate copula covariances from a time series with given levels.
#'  
#' Calculates for each combination of levels \eqn{(\tau_1,\tau_2)}{(tau1,tau2)} 
#' and for all \eqn{k < \code{maxLag}}{k<maxLag} the copula covariances
#' \eqn{Cov(1_{X_0 < \tau_1},1_{X_k < \tau_2})}{Cov(Ind{X0<tau1},Ind{Xk<tau2})}
#' and writes it to \code{values[k]} from its superclass \code{\link{LagOperator}}.
#' 
#' For each lag \code{k = 0, ..., maxLag} and combination of levels
#' \eqn{(\tau_1, \tau_2)}{(tau1, tau2)} from \code{levels.1 x levels.2} the
#' statistic
#' \deqn{\frac{1}{n} \sum_{t=1}^{n-k} ( I\{\hat F_n(Y_t) \leq \tau_1\} - \tau_1) ( I\{\hat F_n(Y_{t+k}) \leq \tau_2\} - \tau_2)}
#' is determined and stored to the array \code{values}.
#' 
#' Currently, the implementation of this class allows only for the analysis of
#' univariate time series.
#' 
#' @name ClippedCov-class
#' @aliases ClippedCov
#'
#' @keywords S4-classes
#' 
################################################################################

setClass(
  Class = "ClippedCov",
  contains = "LagOperator"
)

#' @importFrom stats quantile
#' @importFrom stats acf
setMethod( 
  f = "initialize",
  signature = "ClippedCov",
  definition = function(.Object, Y, maxLag, levels.1, levels.2, isRankBased, positions.boot, B) {
    
    .Object@maxLag = maxLag
    .Object@levels.1 = levels.1
    .Object@levels.2 = levels.2
    .Object@Y = Y
    .Object@positions.boot <- positions.boot
    .Object@B <- B
    
    n <- length(Y)
    ln.1 <- length(levels.1)
    ln.2 <- length(levels.2)
    levels.all = union(levels.1,levels.2)
    ln = length(levels.all)
    
    if(isRankBased){
      if(((max(levels.all) > 1) | (min(levels.all)<0)))
      {stop("all levels must be in [0,1] for a Ranked based estimation")}
      Q = quantile(Y, probs = levels.all)
    }
    
    
    Clipped <- array(0, dim = c(n, ln, B+1))
    for (l in 1:ln) {
      Clipped[,l,1] <- (Y <= Q[l]) - levels.all[l] 
    }
    if (B > 0) {
      pos.boot <- getPositions(.Object@positions.boot,B)
      for (b in 1:B) {
        Clipped[,,b+1] <- Clipped[pos.boot[,b],,1]
      }
    }
    
    pos.1 = match(levels.1,levels.all)
    pos.2 = match(levels.2,levels.all)
    
    val <- array(dim = c(maxLag+1, max(ln.1,1), max(ln.2,1), B+1))
    
    for (b in 0:B) {
      val[,,,b+1] = array(acf(Clipped[,,b+1], type="covariance", lag.max = maxLag, plot = FALSE, demean = FALSE)$acf[1:(maxLag+1),pos.1,pos.2],
          dim = c(maxLag+1, max(ln.1,1), max(ln.2,1)))
    }
    val <- aperm(val, c(1,3,2,4))
    
    .Object@values = val
    
    return(.Object) 
  })

################################################################################
#' Create an instance of the \code{\link{ClippedCov}} class.
#'
#' @name ClippedCov-constructor
#' @aliases clippedCov
#' @export
#'
#' @keywords Constructors
#'
#' @param Y Time series to calculate the copula covariance from
#' @param maxLag maximum lag between observations that should be used
#' @param levels.1 a vector of numerics that determines the level of clipping
#' @param levels.2 a vector of numerics that determines the level of clipping
#' @param isRankBased If true the time series is first transformed to pseudo data;
#' 										currently only rank-based estimation is possible.
#' @param B number of bootstrap replications
#' @param l (expected) length of blocks
#' @param type.boot A flag to choose a method for the block bootstrap; currently
#'                  two options are implemented: \code{"none"} and \code{"mbb"}
#'                  which means to do a moving blocks  bootstrap with \code{B}
#'                  and \code{l} as specified.
#' 
#' @return Returns an instance of \code{ClippedCov}.
#'
#' @seealso \code{\link{LagOperator}}
#'
#' @examples
#' ccf <- clippedCov(rnorm(200), maxLag = 25, levels.1 = c(0.1,0.5,0.9))
#' dim(getValues(ccf))
#' #print values for levels (.5,.5)
#' plot(ccf, maxLag = 20)

################################################################################
clippedCov <- function( Y,
    maxLag = length(Y) - 1,
    levels.1 = c(.5),
    levels.2 = levels.1,
    isRankBased = TRUE,
    B = 0,
    l = 0,
    type.boot = c("none","mbb")){
  
  if (!isRankBased) {
    stop("non rank-based version currently not available")
  }

  if(!(maxLag < length(Y)))
  {maxLag = length(Y) - 1
   warning("maxLag must be smaller then length of dataset, set to maximum")}
  
  if (!((is.vector(levels.1) && is.numeric(levels.2))&&(is.vector(levels.1) && is.numeric(levels.2)))) {
    stop("'levels' needs to be specified as a vector of real numbers")
  }
  
  type.boot <- match.arg(type.boot, c("none","mbb"))[1]
  switch(type.boot,
      "none" = {
        bootPos <- movingBlocks(length(Y),length(Y))},
      "mbb" = {
        bootPos <- movingBlocks(l,length(Y))}
  )
  
  obj = new(
      Class = "ClippedCov",
      Y = Y,
      maxLag = maxLag,
      levels.1 = levels.1,
      levels.2 = levels.2,
      isRankBased = isRankBased,
      B = B,
      positions.boot = bootPos
  )
      
  return(obj)
}