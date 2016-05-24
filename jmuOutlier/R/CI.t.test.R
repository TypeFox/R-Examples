CI.t.test <-
function(x, conf.level=0.95, fpc=1) {
  # Generation of a one-sample two-sided confidence interval on the population mean 
  #    based on the t-test, when sampling without replacement.
  # 'x':  Numeric vector of observations.
  # 'conf.level':  Confidence level of the interval.
  # 'fpc' is the finite population correction, and is used when sampling without replacement.
  # Note: \code{fcp} is typically defined as \code{1-n/N}, 
  #       where \code{n} is the sample size, and \code{N} is the population size.
  # example:  Sample 43 observations from a list of 200 numbers, and compute the 95% confidence interval.
  #           pop = sqrt(1:200) ; x1 = sample( pop, 43 ) ; list(sort(x1))
  #           fpc = 1 - length(x1) / length(pop) ; CI.t.test( x1, fpc=fpc )
  # example:  Sample 14 observations from a Normal(mean=50, sd=5) distribution,
  #              and compute the 90% confidence interval.
  #           x2 = sample( 14, 50, 5 ) ; list(sort(x2)) ; CI.t.test( x2, 0.9 )
  if (!is.numeric(x))  stop("'x' must be numeric.")
  if (length(x)<2)  stop("'x' must contain at least two numbers.")
  if (!is.numeric(conf.level) | length(conf.level)!=1)  stop("'conf.level' must be a number between 0 and 1.")
  if (conf.level <=0 | conf.level >= 1)  stop("'conf.level' must be a number between 0 and 1.")
  if (!is.numeric(fpc))  stop("'fpc' must be numeric.")
  if (fpc <=0 | fpc > 1)  stop("'fpc' must be between 0 and 1.")
  if (fpc==1) {
    y <- as.vector(t.test(x, conf.level=conf.level)$conf.int) }
  else  {
    y <- rep(NA, 2)
    y[1] <- mean(x) + qt( ((1-conf.level)/2), (length(x)-1) ) *
                sqrt( var(x)*fpc/length(x) ) 
    y[2] <- mean(x) - qt( ((1-conf.level)/2), (length(x)-1) ) *
                sqrt( var(x)*fpc/length(x) ) 
  }
  return(y)
}
