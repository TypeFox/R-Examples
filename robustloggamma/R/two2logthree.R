two2logthree <- function(x) {
  if (is.null(x$mu))
    x$mu <- x$shape*x$scale
  c <- 1
  ans <- list()
  ans$mu  <- log(x$scale)+log(x$shape)/c
  ans$sigma <- ans$lambda <- 1/sqrt(x$shape)
  ans$eta <- Exp.response(ans$mu,ans$sigma,ans$lambda)
 
  ans$muse  <- sqrt(avarmu(x$shape)) # you have to divide it by sqrt(n)
  ans$sigmase <- ans$lambdase <- sqrt(avarlambda(x$shape)) # you have to divide it by sqrt(n)
  ans$etase <- sqrt(avareta(x$sigma,x$shape)) # you have to divide it by sqrt(n)
  return(ans)
}
