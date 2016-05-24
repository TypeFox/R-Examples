##### finds a single root to the truncation region function somewhere to the right of 'start'
##### input:
#####     - q, r, s = Coefficients of the truncation region function
#####     - start = Leftmost starting point
#####     - verbose = Print output along the way
find.root <-
function(q, r, s, start, verbose=FALSE){
  right.lim = find.where.sign.changes(q, r, s, start=start, verbose=verbose) # will find a place where we go positive
  if(verbose){print(paste('Found change of sign at: ', right.lim, '; sign: ', sign(truncation.region.function(right.lim, q, r, s)), sep=''), quote=FALSE)}
  x.star = uniroot (f=function(x){truncation.region.function(x, q, r, s)}, lower=start, upper=right.lim)$root
  if(verbose){print(paste('Root at: ', x.star, sep=''), quote=FALSE)}
  x.star
}
