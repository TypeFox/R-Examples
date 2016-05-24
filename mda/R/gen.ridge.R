gen.ridge <-
  function (x, y, weights, lambda = 1, omega, df, ...) 
{
  if (missing(df) && lambda <= .Machine$double.eps) 
    return(polyreg(x, y))
  d <- dim(x)
  mm <- apply(x, 2, mean)
  x <- scale(x, mm, FALSE)
  simple <- if (missing(omega)) 
    TRUE
  else FALSE
  if (!simple) {
    if (!all(match(c("values", "vectors"), names(omega), 
                   FALSE))) 
      stop("You must supply an  eigen-decomposed version of omega")
    vals <- pmax(sqrt(.Machine$double.eps), sqrt(omega$values))
    basis <- scale(omega$vectors, FALSE, vals)
    x <- x %*% basis
  }
  svd.x <- svd(x)
  dd <- svd.x$d
  if (!missing(df))lambda=df.gold(dd,df)
  df=sum(dd^2/(dd^2 + lambda))
  y <- ( t(t(y) %*% svd.x$u) * dd)/(dd^2 + lambda)
  coef <- svd.x$v %*% y
  fitted <- x %*% coef
  if (!simple) 
    coef <- basis %*% coef
  structure(list(fitted.values = fitted, coefficients = coef, 
                 df = df, lambda = lambda, xmeans = mm), class = "gen.ridge")
}
