createWeights2 <- function(I,H){
   
  M <- H[1]
  H <- H[2]
  
  # I number of input neurons
  # M number of hidden neurons in first hidden layer
  # H number of hidden neurons in second hidden layer
  
  # return random weights
  return(new("Weights2", alpha = rnorm(1), alpha_1m = rnorm(M), alpha_2h = rnorm(H), w_h = rnorm(H), 
             q_mh = matrix(nrow=M,ncol=H, data=rnorm(M*H)),w_im = matrix(nrow=I,ncol=M, data=rnorm(I*M))))
  
}