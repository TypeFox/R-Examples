# p = c + (1 - c)ilogit[exp(a*(thet - b))]

p.brm <-
function(theta, params){
  
# Call the appropriate c-function:
  p <- .Call("pbrm", theta, params)
  
# Matrix if n_ppl > 1, vector if n_ppl = 1:
  return(p)

} # END p.brm FUNCTION