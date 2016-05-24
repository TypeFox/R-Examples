residuals.Mort2Dsmooth <-
function(object, 
                                   type = c("deviance",
                                     "pearson",
                                     "anscombe",
                                     "working"),
                                   ...){
  type <- match.arg(type)
  r <- object$residuals
  Z <- object$Z
  fitted.values <- object$fitted.values
  w <- object$w
  res <- switch(type,
                deviance = r, 
                pearson  = (Z - fitted.values) / sqrt(fitted.values),
                anscombe = (3/2) * (Z^(2/3) - fitted.values^(2/3)) / fitted.values^(1/6),
                working  = (Z - fitted.values)
                )
  res <- matrix(res,
                object$m,
                object$n,
                dimnames=list(object$x, object$y))
  res
}
