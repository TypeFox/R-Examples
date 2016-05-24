#' Summary statistics for class ELT.
#'
#' Summary statistics for class ELT.
#'
#' @param object An object of class \code{ELT}. Data frame containing two numeric columns. The column \code{Loss} contains the expected losses from each single occurrence of event. The column \code{Rate} contains the arrival rates of a single occurrence of event. 
#' @param t Scalar representing the time period of interest. The default value is \code{t} = 1.
#' @param theta Scalar containing information about the variance of the Gamma distribution: \eqn{sd[X] = x * }\code{theta}. The default value is \code{theta} = 0: the loss associated to an event is considered as a constant.
#' @param cap Scalar representing the financial cap on losses for a single event, i.e. the maximum possible loss caused by a single event. The default value is \code{cap} = Inf.
#' @param ... additional arguments affecting the summary produced.

#' @return A list containing the data summary, and the means and the standard deviations of \eqn{N}, \eqn{Y}, and \eqn{S}. 

#' @method summary ELT

#' @export
#' @examples
#' data(UShurricane)
#' summary(ELT(UShurricane))

summary.ELT<- function(object, theta = 0, cap = Inf, t = 1, ...){
	
		stopifnot(inherits(object, "ELT"),
             is.scalar(t), t > 0,
             is.scalar(cap), cap > 0,
             is.scalar(theta), theta >= 0)
  
  mom <- getMoments(object$Loss, theta = theta, cap = cap, maxmom = 2)
  lambda <- sum(object$Rate)
  EY <- (crossprod(object$Rate, mom)) / lambda
  ES <- lambda * t * EY
  
  res <- c(lambda * t, sqrt(lambda * t), EY[1], sqrt(EY[2] - EY[1]^2), ES[1], sqrt(ES[2]))
  
  names(res) <- c(paste("E[N(", t, ")]", sep = ""), paste("Sd[N(", t, ")]", sep = ""), "E[Y]", "Sd[Y]", paste("E[S(", t, ")]", sep = ""), paste("Sd[S(", t, ")]", sep = ""))
  
  print(list("Data Summary" = summary.data.frame(object, ...), "Summary Statistics" = res))
  
}
