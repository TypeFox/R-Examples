##### returns the negative for a function known to be decreasing everywhere to the right of 'start'
##### input:
#####     - q, r, s = Coefficients of the truncation region function
#####     - start = Leftmost point at which we know the derivative is negative
#####     - verbose = Print output along the way
interval.if.derivative.always.negative <-
function(q, r, s, start, verbose=FALSE){
  f.start = truncation.region.function(start, q, r, s)
  if (f.start < 0){ # always negative - so return full half line
    if (verbose){print ('Starts negative; keeps decreasing. Returning half line.', quote=FALSE)}
    return (Intervals(c(0, Inf)))
  }else{ # eventually turns negative
    x.star = find.root (q, r, s, start=start, verbose=verbose)
    return (Intervals(c(x.star, Inf)))
  }
}
