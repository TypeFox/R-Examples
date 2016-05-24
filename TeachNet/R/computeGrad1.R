computeGrad1 <- function(x,y,I,H,weights,f,f_d,m_f){
  # computes gradient for one observation pair
  # x are the covariates
  # y the class
  # I number of input variable
  # H number of hidden neuron in first hidden layer
  # weights for this step
  # f the acctivation function
  # f_d derivate of acctivation function
  # m_f function to calculate m
  
  # calculate interim vector r
  r <- c(1:H)
  r <- vapply(r, function(h) weights@alpha_h[h] + sum(weights@w_ih[ ,h]*x),1)
  
  # calculate interim value z
  z <- weights@alpha + sum(weights@w_h*f(r))
  
  # calculate interim value 2*[y^ -y]
  m <- m_f(z,y)
  
  # create gradient
  grad <- new("Weights", alpha = 0, alpha_h = c(1:H), w_h = c(1:H), w_ih = matrix(nrow=I,ncol=H, data=0))
  
  # Part 1
  grad@alpha <- m*f_d(z)

  
  # Part 2 and 3
  grad@w_h <-sapply(grad@w_h, function(h) grad@alpha*f(r[h]))
  grad@alpha_h <- sapply(grad@alpha_h, function(h) grad@alpha*weights@w_h[h]*f_d(r[h]))
  
  # Part 4
  grad@w_ih <- outer(1:I, 1:H , FUN=function(i,h) grad@alpha_h[h]*x[i] )
  
  return(grad)
  
}# end function