createWeights1 <- function(I,H){
  # I number of input neurons
  # H number of hidden neurons
  
  # return random weights
  return(new("Weights", alpha = rnorm(1), alpha_h = rnorm(H), w_h = rnorm(H), 
            w_ih = matrix(nrow=I,ncol=H, data=rnorm(I*H))))
  
}