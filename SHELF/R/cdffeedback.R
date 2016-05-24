#' Feedback for the elicited distribution of the population CDF
#' 
#' Report the median and 100(1-alpha)\% credible interval for point on the population CDF
#' 
#' Denote the uncertain population CDF by \deqn{P(X \le x | \mu, \sigma^2),}where \eqn{\mu}
#' is the uncertain population median and \eqn{\sigma^(-2)} is the uncertain population precision.
#' Feedback can be reported in the form of the median and 100(1-alpha)\% credible interval for
#' (a) an uncertain probability \eqn{P(X \le x | \mu, \sigma^2)}, where \eqn{x} is a specified 
#' population value and (b) an uncertain quantile \eqn{x_q} defined by \eqn{P(X \le x_q | \mu, \sigma^2) = q}, where \eqn{q} is a specified 
#' population probability.  
#'
#' @param medianfit The output of a \link{fitdist} command following elicitation
#'  of the expert's beliefs about the population median.
#' @param precisionfit The output of a \link{fitprecision} command following elicitation
#'  of the expert's beliefs about the population precision.
#' @param quantiles A vector of quantiles \eqn{q_1, \ldots,q_n} required for feedback
#' @param vals A vector of population values \eqn{x_1,\ldots,x_n} required for feedback
#' @param alpha The size of the 100(1-alpha)\% credible interval
#' @param median.dist The fitted distribution for the population median. Can be one of \code{"normal"},
#'  \code{"lognormal"} or \code{"best"}, where \code{"best"} will select the best fitting out of 
#'  normal and lognormal.
#' @param precision.dist The fitted distribution for the population precision. Can either be \code{"gamma"}
#'  or \code{"lognormal"}. 
#' @param n.rep The number of randomly sampled CDFs used to estimated the median
#'  and credible interval.
#'  
#' @return Fitted median and 100(1-alpha)\% credible interval for population 
#'  quantiles and probabilities.
#'  
#'  \item{$quantiles}{Each row gives the fitted median 
#'  and 100(1-alpha)\% credible interval for each uncertain population quantile 
#'  specified in \code{quantiles}: the fitted median 
#'  and 100(1-alpha)\% credible interval for the value of \eqn{x_{q_i}} where 
#'  \eqn{P(X\le x_{q_i} | \mu, \sigma^2) = q_i.}} 
#'  \item{$probs}{Each row gives the fitted median 
#'  and 100(1-alpha)\% credible interval for each uncertain population probability 
#'  specified in \code{probs}: the fitted median 
#'  and 100(1-alpha)\% credible interval for the value of
#'  \eqn{P(X\le x_i | \mu, \sigma^2).} }
#'
#' @examples 
#' \dontrun{
#' prfit <- fitprecision(interval = c(60, 70), propvals = c(0.2, 0.4), trans = "log")
#' medianfit <- fitdist(vals = c(50, 60, 70), probs = c(0.05, 0.5,  0.95), lower = 0)
#' cdffeedback(medianfit, prfit, quantiles = c(0.01, 0.99),
#'             vals = c(65, 75), alpha = 0.05, n.rep = 10000)
#'  }  
#' @export

cdffeedback <- function(medianfit, precisionfit, quantiles = c(0.05, 0.95), 
                        vals = NA, alpha = 0.05, median.dist = "best", 
                        precision.dist = "gamma", n.rep = 10000){
  
  if(precision.dist!="gamma" & precision.dist!="lognormal"){
    stop('precision.dist must equal one of "gamma" or "lognormal"')
  }
  
  f <- getdists(precisionfit$transform)
  mediandist <- getmediandist(medianfit, median.dist)
  
  musample <- mediandist$rand(n.rep, mediandist$m, mediandist$s)
  mumatrix <- matrix(musample, n.rep, length(quantiles))
  
  if(precision.dist == "gamma"){
    sigmasample <- sqrt(1 / rgamma(n.rep, precisionfit$Gamma[[1]], 
                                   precisionfit$Gamma[[2]]))
  }
  
  if(precision.dist == "lognormal"){
    sigmasample <- sqrt(1 / rlnorm(n.rep, precisionfit$Log.normal[[1]], 
                                   precisionfit$Log.normal[[2]]))
  }
  
  sigmamatrix <- matrix(sigmasample, n.rep, length(quantiles)) 
  
  quantilematrix <- matrix(quantiles, n.rep, length(quantiles), byrow = T)
  quantilesample <- f$quan(quantilematrix, f$trans(mumatrix), sigmamatrix) 
  quantilesample <- apply(quantilesample, 2, sort)
  
  index <- c(ceiling(alpha / 2 * n.rep), round(0.5 * n.rep),  floor((1 - alpha / 2)*n.rep))
  quantileinterval <- quantilesample[index, ]
  
  rownames(quantileinterval) <- c(alpha / 2, 0.5, 1 - alpha/2)
  colnames(quantileinterval) <- quantiles
  
  if(!is.na(vals[1])){
    mumatrix <- matrix(musample, n.rep, length(vals))
    sigmamatrix <- matrix(sigmasample, n.rep, length(vals))
    valmatrix <- matrix(vals, n.rep, length(vals), byrow = T)
    probsample <- f$cdf(valmatrix, f$trans(mumatrix), sigmamatrix) 
    probsample <- apply(probsample, 2, sort)
    probinterval <- probsample[index, ]
    rownames(probinterval) <- c(alpha / 2, 0.5, 1 - alpha/2)
    colnames(probinterval) <- vals
    return(list(quantiles = t(quantileinterval), probs = t(probinterval)))}else{
      return(list(quantiles = t(quantileinterval)))
    }
}
