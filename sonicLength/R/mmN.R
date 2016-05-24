mmN <-
function(probs,
         max.N=log(.Machine$double.eps)/log(max(pmin(1-probs,1-1e-7))),
         tol=sqrt(.Machine$double.eps))
{

  ## Purpose: Create a function to map number of lengths to estimate
  ##            number of IBD sonicants

  ## probably mask this from user view - its a utility function that
  ## should be called by some more vectorized version.  it is slow,
  ## but only need to call it for unique values, perhaps like this:
  ##    res <- mmN( pr.len )
  ##    x2 <- ave(x,x,FUN=function(m) res(m[1]))
  
  ## ----------------------------------------------------------------------
  ## Arguments: probs - probabilities of various lengths
  ##            max.N - largest N.val to use 
  ##            tol - objective needs to be no more than this
  ## ----------------------------------------------------------------------
  ## Author: Charles Berry, Date:Sun Mar 20 2011

  f <- function(x,probs,m) (sum(1 - (1-probs)^x ) - m)
  function(m) {
    stopifnot ( length(m)==1 )
    stopifnot ( m < sum( probs > 0 ))
    ## numerically f(m,probs,m) may be > 0
    tryval <- f(m,probs,m)
    stopifnot( tryval< tol )
    res <-
      if (tryval>0) m else uniroot(f,c(m,max.N), m=m, probs=probs)$root
    
    res
  }
}
