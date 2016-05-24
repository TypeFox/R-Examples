#' @importFrom stats qnorm
getCI <- function(ts, var, fisher, ci=0.95){
  if(!fisher){
    lower <- ts - abs(qnorm(0.5*(1-ci)))*sqrt(var)
    upper <- ts + abs(qnorm(0.5*(1-ci)))*sqrt(var)  
  } else {
    ts_f <- log( (1+ts)/(1-ts) )
    var_f <- var*(2/(1-ts^2))^2
    lower_f <- ts_f - abs(qnorm(0.5*(1-ci)))*sqrt(var_f)
    upper_f <- ts_f + abs(qnorm(0.5*(1-ci)))*sqrt(var_f)
    lower <- (exp(lower_f)-1)/(1+exp(lower_f))
    upper <- (exp(upper_f)-1)/(1+exp(upper_f))
  }
  return(c(lower, upper))
}

