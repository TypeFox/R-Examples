#### evaluates the function f(x) = q.sqrt(x) + r.sqrt(1+x) + s at start - stores sign
#### then updates by multiplying and adding until we get to opposite sign
#### input:
####    - q, r, s = Truncation region coefficients
####    - start = Where do we start (left)?
####    - add, multiply = Constants defining the way in which we update the next guess in our search for something with opposite sign
####    - verbose = Print output along the way
find.where.sign.changes <-
function(q, r, s, start, add=1, multiply=2, verbose=FALSE){
  initial.sign = sign(truncation.region.function(start, q, r, s))
  if(verbose){print (paste('Initial sign at ', start, ': ', initial.sign, sep=''), quote=FALSE)}
  now.try = start
  
  while (TRUE){
    now.try = multiply*now.try + add
    if (verbose){print (paste('Now trying: ', now.try, sep=''), quote=FALSE)}
    new.sign = sign(truncation.region.function(now.try, q, r, s))
    
    if (new.sign != initial.sign){break}
  }
  
  return (now.try)
}
