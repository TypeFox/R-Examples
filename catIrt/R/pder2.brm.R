# p'' = a*[1 - exp(a*(thet - b))]*q^(1)*pder1
# - p^(1): 2PL probability
# - q^(1): 1 - p^(1)

pder2.brm <-
function(theta, params){
  
# Call the appropriate c-function:
  pder2 <- .Call("pder2brm", theta, params)
  
# Matrix if n_ppl > 1, vector if n_ppl = 1:
  return(pder2)
     
} # END pder2.brm FUNCTION

