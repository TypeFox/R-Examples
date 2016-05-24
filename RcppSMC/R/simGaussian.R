simGaussian <- function(len = 250)
{
  sim <- list()
  sim$state <- cumsum(rnorm(len))
  sim$data  <- sim$state + rnorm(len)

  invisible(sim)
}
