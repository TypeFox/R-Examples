# p' = a*(p.ku*q.ku - p.kl*q.kl)

pder1.grm <-
function(theta, params){
  
# Call the appropriate c-function:
  pder1 <- .Call("pder1grm", theta, params)
  
# Matrix of (n_ppl * n_cat) x n_it:
  return(pder1)
 
} # END pder1.grm FUNCTION

