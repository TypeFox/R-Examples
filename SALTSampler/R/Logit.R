Logit <- function(p, logp = FALSE) {
  #Return logit of probabilities p. If logp == TRUE, the first argument 
  #contains natural logs of probabilites which helps preserve accuracy
  #for probailities near 1
  
  #Check Input
  if (!exists("p")) {
    stop("p is not defined")
  }
  
  if (!is.numeric(p)) {
    stop("p is not numeric")
  }
  if(logp == FALSE & (any(p > 1) || any(p < 0))) {
    stop("p is not a probability")
  }
  if(logp == TRUE & (any(p > log(1)) || any(p < log(0)))) {
    stop("p is not a probability")
  }
  
  #Two different ways to write logit(p) on natural scale
  #Which method is used depends on which is more numerically 
  #accurate
  out <- numeric(length(p))
  if (!logp) {
    ok <- (p < 1/2)
    if (any(ok)) out[ok] <- log(p[ok]) - log1p(-p[ok])
    if (any(!ok)) out[!ok] <- log(p[!ok]/(1 - p[!ok]))
    return(out)
  }
  
  #Two different ways to write logit(p) on log scale
  #Which method is used depends on which is more numerically 
  #accurate
  lp <- p 
  ok <- (lp < log(1/2))
  if (any(ok)) out[ok] <- lp[ok] - log1p(-exp(lp[ok]))
  if (any(!ok)) out[!ok] <- lp[!ok] - log(-expm1(lp[!ok]))
  
  #Return logit values
  return(out)
}
