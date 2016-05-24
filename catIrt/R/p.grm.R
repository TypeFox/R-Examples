# p = ilogit[exp(a*(thet - b.ku))] - ilogit[exp(a*(thet - b.kl))]

p.grm <-
function(theta, params){
  
# Call the appropriate c-function:
  p <- .Call("pgrm", theta, params)
  
# Matrix of (n_ppl * n_cat) x n_it:
  return(p)

} # END p.grm FUNCTION
