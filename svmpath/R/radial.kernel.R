"radial.kernel" <-
  function(x, y=x, param.kernel = 1/p,...)
{

  ###Note param.kernel is now gamma 
  n <- nrow(x)
  m <- nrow(y)
  p <- ncol(x)
  normx <- drop((x^2) %*% rep(1, p))
  normy <- drop((y^2) %*% rep(1, p))
  a <- x %*% t(y)
  a <- (-2 * a + normx) + outer(rep(1, n), normy, "*")
  exp( - a* param.kernel)
}
