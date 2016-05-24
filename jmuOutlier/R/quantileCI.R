quantileCI <-
function( x, probs = 0.5, conf.level = 0.95 ) { 
  # Produces exact confidence intervals on quantiles corresponding to the given probabilities,
  #    based on the binomial test. 
  # 'x': Numeric vector of observations.
  # 'probs': Numeric vector of cumulative probabilities.
  # 'conf.level':  Confidence level of the interval.
  # Note:  If \code{probs} is 0.5 (default), then a confidence interval is constructed on the population median.
  # Example:  # Sample 20 observations from an Exponential distribution with mean=10.
  #           # Construct 90% confidence intervals on the 25th, 50th, and 75th percentiles.  
  #           x = rexp( 20, 0.1 ) ;   list(sort(x)) ;   quantileCI( x, c( 0.25, 0.5, 0.75 ), 0.9 )
  if (!is.numeric(x))  stop("'x' must be numeric.")
  if (length(x)<2)  stop("'x' must contain at least two numbers.")
  if (!is.numeric(probs))  stop("'probs' must be numeric.")
  if (length(x)<2)  stop("'x' must contain at least two numbers.")
  if (!is.numeric(probs))  stop("'probs' must be numeric.")
  if (min(probs)<=0 | max(probs)>=1)  stop("'probs' must be between 0 and 1.")
  if (!is.numeric(conf.level) | length(conf.level)!=1)  stop("'conf.level' must be a number between 0 and 1.")
  if (conf.level <=0 | conf.level >= 1)  stop("'conf.level' must be a number between 0 and 1.")
  CI = NULL
  for ( k in 1:length(probs) )   {
     delta <- (max(x)-min(x))/1e10
     xgrid <- c(x, x+delta, x-delta)
     value.in.CI <- rep(NA, length(xgrid))
     for (iii in 1:length(xgrid)) {
        x1 <- c( sum( x<xgrid[iii] ), sum( x>xgrid[iii] ) ); n <- sum(x1)
        value.in.CI[iii] <-
           binom.test(x1, n, probs[k], alternative = "two.sided", conf.level)$p.value>=1-conf.level
     }
     CI <- rbind( CI,  c( min( xgrid[value.in.CI] ), max( xgrid[value.in.CI] ) ) )
  }
  colnames( CI ) = c("lower","upper") ;  rownames( CI ) = probs 
  return( CI )
}
