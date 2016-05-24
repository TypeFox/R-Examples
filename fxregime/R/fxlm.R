## Frankel-Wei linear regression (simple wrapper for lm())
fxlm <- function(formula, data, ...) {
  cl <- match.call()
  if(missing(formula)) formula <- colnames(data)[1]
  if(is.character(formula)) formula <- as.formula(paste(formula,
    "~", paste(colnames(data)[colnames(data) != formula], collapse = " + ")))
  rval <- lm(formula, data = as.data.frame(data), ...)
  rval$call <- cl
  rval$index <- index(data)
  class(rval) <- c("fxlm", "lm")
  return(rval)
}

## and corresponding coefficients and estimating functions
coef.fxlm <- function(object, ...) {
  rval <- NextMethod(object)
  rval <- c(rval, "(Variance)" = mean(residuals(object)^2))
  return(rval)
}

## estfun.fxlm <- function(x, ...) {
##   res <- residuals(x)
##   ef <- NextMethod(x)
##   sigma2 <- mean(res^2)
##   rval <- cbind(ef, (res^2 - sigma2))
##   colnames(rval) <- c(colnames(ef), "(Variance)")
##   if(!inherits(rval, "zoo"))
##     rval <- zoo(rval, index(x))
##   return(rval)
## }
## 
## bread.fxlm <- function(x, ...) {
##   br <- NextMethod(x)
##   sigma2 <- coef(x)["(Variance)"]
##   br <- rbind(cbind(br, "(Variance)" = 0), "(Variance)" = 0)
##   br[nrow(br), ncol(br)] <- 1
##   br
## }

estfun.fxlm <- function(x, ...) {
  res <- residuals(x)
  ef <- NextMethod(x)
  sigma2 <- mean(res^2)
  rval <- cbind(ef/sigma2, (res^2 - sigma2)/(2 * sigma2^2))
  colnames(rval) <- c(colnames(ef), "(Variance)")
  if(!inherits(rval, "zoo"))
    rval <- zoo(rval, index(x))
  return(rval)
}

bread.fxlm <- function(x, ...) {
  br <- NextMethod(x)
  sigma2 <- coef(x)["(Variance)"]
  ef <- estfun(x)
  b1 <- solve(br)/sigma2
  b2 <- colMeans(ef[,-ncol(ef)]/sigma2)
  b3 <- 0.5/(sigma2^2)
    
  br <- rbind(cbind(b1, "(Variance)" = b2), "(Variance)" = c(b2, b3))
  solve(br)
}

## index/time extraction
time.fxlm <- index.fxlm <- function(x, ...) x$index

## test linear hypotheses
linearHypothesis.fxlm <- function(model, ...) {
  class(model) <- class(model)[-1]
  NextMethod()
}
