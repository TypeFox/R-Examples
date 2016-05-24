residuals.Mort1Dsmooth <-
function(object,
                                   type = c("deviance",
                                     "pearson",
                                     "anscombe",
                                     "working"), ...){
  type <- match.arg(type)
  r <- object$residuals
  y <- object$y
  fitted.values <- object$fitted.values
  w <- object$w
  res <- switch(type,
                deviance = r, 
                pearson  = (y - fitted.values) / sqrt(fitted.values),
                anscombe = (3/2) * (y^(2/3) - fitted.values^(2/3)) / fitted.values^(1/6),
                working  = (y - fitted.values)
                )
  res
}
