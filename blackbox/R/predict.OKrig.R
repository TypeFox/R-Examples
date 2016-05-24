predict.OKrig <- function (object, x = NULL,...)  # ... for consistency with generic
{
  if (is.null(x)) { ## call 'predict(fitobject)'
    x <- object$x
  } else if ( class(x)=="numeric" ) {
    x <- t(x) ## row vector in matrix class
  } else if ( ! is.matrix(x) ) { ## 2D object of the wrong class
    x <- t(t(x))
  }
  if (nrow(x)==1) { ## single point
    temp <- Cpredict(x, object)
  } else {
    temp <- apply(x, 1, function(v) {Cpredict(v, object)})
  }
  return(c(temp))
}
