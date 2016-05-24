# robust estimator of standard deviation for Gaussian data with jumps
"sdrobnorm" <-
function (x, p = c(0.25, 0.75), lag = 1)
{
  as.numeric( diff(quantile( diff(x, lag = lag), p)) / diff(qnorm(p)) / sqrt(2) )
}
