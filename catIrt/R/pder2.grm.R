# p'' = a*[1 - exp(a*(thet - b.ku))]*q.ku*pder1.ku - a*[1 - exp(a*(thet - b.kl))]*q.kl*pder1.kl

pder2.grm <-
function(theta, params){
  
# Call the appropriate c-function:
  pder2 <- .Call("pder2grm", theta, params)
  
# Matrix of (n_ppl * n_cat) x n_it:
  return(pder2)
     
} # END pder2.grm FUNCTION

