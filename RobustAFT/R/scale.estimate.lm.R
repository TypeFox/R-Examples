scale.estimate.lm <-
function(object,...)
 {
  rdf <- object$df.resid
  resid <- object$residuals
  wt <- object$weights 
  n <- length(resid)
  p <- object$rank
  if (is.null(rdf)) rdf <- n-p
  if(!is.null(wt)) {
                wt <- wt^0.5
                resid <- resid * wt
                excl <- wt == 0
                rdf <- rdf - sum(excl)}
  if (n > p) z <- sqrt(sum(resid^2)/rdf)
  else z <- NA  
  class(z) <- "lm.i"
 z}

