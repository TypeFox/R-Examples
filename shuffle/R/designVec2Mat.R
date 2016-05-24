designVec2Mat <-
function(desvec){ 
  #########
  # Convert from a design vector (a map t->j)
  # into a design matrix (X_t,j=1 is t->j)
  # if desvec[t]=0 they are ignored
  #
  # If (a_1...a_m) is a vector of effects,
  # then Y_t = desmat_t*a + eps_t and Y_t = desvec(t) + eps_t 
  #########
  desmat = matrix(0,nrow=length(desvec),ncol = max(desvec))
  for (j in 1:max(desvec)){
    desmat[,j] = (desvec==j)
  }
  return(desmat)
}
