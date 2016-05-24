# p' = (1 - c)*a*p^(1)*q^(1)
# - p^(1): 2PL probability.
# - q^(1): 1 - p^(1).

pder1.brm <-
function(theta, params){
  
# Call the appropriate c-function:
  pder1 <- .Call("pder1brm", theta, params)
  
# Matrix if n_ppl > 1, vector if n_ppl = 1:
  return(pder1)
 
} # END pder1.brm FUNCTION
