#' Fit a distribution to judgements about a population precision
#' 
#' Takes elicited probabilities about proportion of a population
#' lying in a specfied interval as inputs, converts the judgements into probability
#' judgements about the population precision, and fits gamma and lognormal distributions
#' to these judgements using the \link{fitdist} function. 
#' 
#' The expert provides a pair of probability judgements
#'  \deqn{P(\theta < \theta_1 ) = p_1,} and \deqn{P(\theta < \theta_2) = p_2,}
#'  where \eqn{\theta} is the proportion of the population that lies in the interval
#'  \eqn{[k_1, k_2]}. The judgements are made conditional on the population median 
#'  equalling \eqn{k_1}. Note that, unlike the \link{fitdist} command, a 'best fitting'
#'  distribution is not reported, as the distributions are fitted to two elicited
#'  probabilities only.    
#' 
#' @param interval A vector specifying the endpoints of an interval \eqn{[k_1, k_2]}.  
#' @param propvals A vector specifying two values \eqn{\theta_1, \theta_2} for the proportion.
#' @param propprobs A vector specifying two probabilities \eqn{p_1, p_2}.
#' @param trans A string variable taking the value \code{"identity"}, \code{"log"} or
#' \code{"logit"} corresponding to whether the population distribution is normal, lognormal
#' or logit-normal respectively.
#' @param pplot Plot the population distributions with median set at \eqn{k_1}
#' and precision fixed at the two elicited quantiles implied by \code{propvals} 
#' and \code{propprobs}.
#' @param fontsize Font size used in the plots.
#' 
#' @return 
#' \item{Gamma}{Parameters of the fitted gamma distribution. Note that E(precision) =
#' shape / rate.} 
#' \item{Log.normal}{Parameters of the fitted log normal
#' distribution: the mean and standard deviation of log precision.}
#' \item{vals}{The elicited values \eqn{\theta_1, \theta_2}}
#' \item{probs}{The elicited probabilities \eqn{p_1, p_2}}
#' \item{limits}{The lower and upper limits specified by each expert (+/- Inf
#' if not specified).}
#' \item{transform}{Transformation used for a normal population distribution.}

#' 
#'    
#' @examples 
#' \dontrun{
#' fitprecision(interval=c(60, 70), propvals=c(0.2, 0.4), trans = "log")
#'   }
#' @export

fitprecision <- function(interval, propvals, 
                         propprobs = c(0.05, 0.95), 
                         trans = "identity", pplot = TRUE,
                         fontsize = 18){
  
  if (max(propvals >=0.5) | min(propvals <=0)){
    stop('propvals must be between 0 and 0.5')
  }
  
  if (trans != "identity" & trans != "log" & trans != "logit"){
    stop('argument trans must be one of "identity", "log" or "logit"')}
  
  if (trans == "identity"){
    precisionvalues <- (qnorm(propvals + 0.5) / (interval[2] - interval[1]))^2
    dens <- dnorm
    quan <- qnorm
    m <- interval[1]
  }
  
  if (trans == "log"){
    precisionvalues <- (qnorm(propvals + 0.5) / (log(interval[2]) - log(interval[1])))^2
    dens <- dlnorm
    quan <- qlnorm
    m <- log(interval[1])
  }
  
  if (trans == "logit"){
    precisionvalues <- (qnorm(propvals + 0.5) / (logit(interval[2]) - logit(interval[1])))^2
    dens <- dlogit
    quan <- qlogit
    m <- logit(interval[1])
  }
  
  precisionfit <- fitdist(vals = precisionvalues, probs = propprobs, lower = 0)
  precisionfit$transform <- trans
  
  if(pplot == TRUE){
    
    s <- sort(1 / sqrt(precisionvalues))
    xl<-quan(0.001, m, s[2])
    xu<-quan(0.999, m, s[2])
    x <- seq(from = xl, to  = xu, length = 200)
    d1 <- dens(x, m, s[2])
    xint <- seq(from = interval[1], to = interval[2], length = 200)
    dint1 <- dens(xint, m, s[2])
    d2 <- dens(x, m, s[1])
    dint2 <- dens(xint, m, s[1])
    df<-data.frame(x=x, d1=d1, d2=d2, xint=xint, dint1=dint1, dint2=dint2)
    theme_set(theme_grey(base_size = fontsize))
    
    pcore <- ggplot(df, aes(x=x, y=d1)) + expand_limits(y = c(0, max(d2))) + labs(y = "")
    p1 <- pcore + geom_line() + 
      geom_area(data = df, aes(x=xint, y = dint1), fill="red", alpha=0.5) +
      labs(title = paste("lower (",propprobs[1],
                         " quantile) proportion = ", propvals[1], sep=""))
    p2 <- pcore + geom_line(aes(x=x, y=d2)) + 
      geom_area(data = df, aes(x=xint, y = dint2), fill="red", alpha=0.5) +
      labs(title = paste("upper (",propprobs[2],
                         " quantile) proportion = ", propvals[2], sep =""))
    multiplot(p1, p2)
  }
  
  
  precisionfit$Normal <- precisionfit$Student.t <- precisionfit$Log.Student.t <- NULL
  precisionfit$best.fitting <- precisionfit$Beta <- precisionfit$ssq <- NULL
  precisionfit$limits <- NULL
  precisionfit$vals <- propvals
  precisionfit
}