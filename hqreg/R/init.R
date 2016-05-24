# Determine vector d used to compute lambda.max
initDeriv <- function(yy, X, method, gamma, c, penalty.factor)
{
  ind <- which(penalty.factor!=0)
  if(length(ind) == ncol(X)) {
    if(method == "huber") {
      d <- ifelse(abs(yy)>gamma, sign(yy), yy/gamma)
    } else if(method == "quantile") {
      d <- ifelse(abs(yy)>gamma, sign(yy), yy/gamma)+c
    } else {
      d <- yy
    }
  } else {
    if(method == "huber") {
      fit <- lm(yy~X[,-ind]+0)
      r <- fit$residuals
      d <- ifelse(abs(r)>gamma, sign(r), r/gamma)
    } else if(method == "quantile") {
      fit <- lm(yy~X[,-ind]+0)
      r <- fit$residuals
      d <- ifelse(abs(r)>gamma, sign(r), r/gamma)+c
    } else {
      fit <- lm(yy~X[,-ind]+0)
      d <- fit$residuals
    } 
  }
  d
}

