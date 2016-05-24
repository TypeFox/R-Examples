# References
#   http://en.wikipedia.org/wiki/Uhlenbeck-Ornstein_process
ou.process <- function(theta, mu=0, sigma=1, initial=mu,
  end=Sys.Date(), start=NULL, obs=NULL)
{
  dates <- getTradingDates(end, start, obs)

  if (is.null(obs)) obs <- length(dates)
  ts <- 1:obs
  series <- initial * exp(-theta *ts) + mu * (1 - exp(-theta * ts))
  noise <- rnorm(obs, sd=sqrt(sigma^2 / (2*theta) * (1-exp(-2*theta*ts)) ) )
  xts(series + noise, order.by=dates)
}

# TODO: Add common interfaces for these processes
