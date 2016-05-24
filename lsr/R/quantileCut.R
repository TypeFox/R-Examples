# file:    quantileCut.R 
# author:  Dan Navarro
# contact: daniel.navarro@adelaide.edu.au
# changed: 13 November 2013

# quantileCut() splits the data x into n equally sized groups. Use with care: 
# it is very frequently the case that people want to break a variable into several
# equally sized groups in order to force a continuous variable into the ANOVA
# framework. Much as I recognise the intuitive appeal of this, I'm not convinced
# it's actually a sensible data analysis strategy *unless* the groups that emerge
# from this automatic-grouping strategy are actually meaningful. I've included the
# function because it's something that many students have requested, but as I say,
# it should be used with care.
quantileCut <- function(x, n, ...) {
  
  if( !is.vector(x) | !is(x,"numeric")) stop( '"x" must be a numeric vector')
  if( length(n) !=1 | !is(n,"numeric")) stop( 'number of bins "n" must be a single number')
  
  p <- seq(0,1,1/n)
  breaks <- quantile( x, p, na.rm=TRUE )
  eps <- (max(x, na.rm=TRUE)-min(x, na.rm=TRUE)) / 1000
  breaks[1] <- breaks[1] - eps
  breaks[n+1] <- breaks[n+1] + eps  
  out <- cut(x, breaks, ...)
  return( out )
}
