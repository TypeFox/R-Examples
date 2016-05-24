simLineart <- function(len = 250)
{
  sim <- list()

  statex <- cumsum(rnorm(len))
  statey <- cumsum(rnorm(len))

  datax  <- statex + rt(len,df=4)
  datay  <- statey + rt(len,df=4)

  sim$state <- matrix(cbind(statex,statey), ncol=2, dimnames=list(NULL, c("x", "y")))
  sim$data  <- matrix(cbind(datax,datay), ncol=2, dimnames=list(NULL, c("x", "y")))

  invisible(sim)
}
