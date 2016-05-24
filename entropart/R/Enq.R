Enq <-
function(n, q)
{
  if (q == 0) {
    return (rep(1, length(n)))
  } else {
    # beta cannot be computed for n-q+1<0 (so warnings must be suppressed) but the value is 0 then
    BetaVersion  <- suppressWarnings(gamma(q)/beta(n-q+1, q))
    BetaVersion[n-q+1<0] <- 0
    # (-1)^n is problematic for long vectors (returns NA for large values). It is replaced by 1-n%%2*2. It is replaced by 1-n%%2*2. (n is rounded if is not an integer)
    ExpNq <- BetaVersion-(1-round(n)%%2*2)*gamma(1+q)*sin(pi*q)/pi/(n+1)   
    return (ExpNq)
  }
}                                
