##### returns the negative for a function known to be increasing everywhere to the right of 'start'
##### input:
#####     - q, r, s = Coefficients of the truncation region function
#####     - start = Leftmost point at which we know the derivative is positive
#####     - verbose = Print output along the way
interval.if.derivative.always.positive <-
function(q, r, s, start, verbose=FALSE){
  f.start = truncation.region.function(start, q, r, s)
  if (f.start > 0){ # always positive - return empty set
    if (verbose){print ('Starts positive; keeps increasing. Returning empty set.', quote=FALSE)}
    return (Intervals(c(0, 0)))
  }else{ # goes positive somewhere
    x.star = find.root (q, r, s, start=start, verbose=verbose)
    return (Intervals(c(0, x.star)))
  }
}
