computeGrad2 <- function(x,y,I,M,H,weights,f,f_d,m_f){
  # computes gradient for one observation pair
  
  # x are the covariates
  # y the class
  # I number of input variable
  # H number of hidden neuron in first hidden layer
  # H number of hidden neuron in second hidden layer
  # weights for this step
  # f the acctivation function
  # f_d derivate of acctivation function
  # m_f function to calculate m
  
  # create gradient
  grad <- new("Weights2", alpha = 0, alpha_1m = c(1:M), alpha_2h = c(1:H), w_h = c(1:H), 
              q_mh = matrix(nrow=M,ncol=H, data=0), w_im = matrix(nrow=I,ncol=M, data=0))
  
  # calculate interim vector r
  r <- c(1:M)
  r <- vapply(r, function(m) weights@alpha_1m[m] + sum(weights@w_im[ ,m]*x),1)
  # calculate interim value z
  z <- c(1:H)
  z <- vapply(z, function(h) weights@alpha_2h[h] + sum(weights@q_mh[ ,h]*f(r)),1)
  
  # calculate interim value s
  s <- weights@alpha + sum(weights@w_h*f(z))
  
  # calculate interim value 2*[y^ -y]
  p <- m_f(s,y)
  
  # Part 1
  grad@alpha <- p*f_d(s)

  # Part 2 and 3
  grad@w_h <-sapply(grad@w_h, function(h) grad@alpha*f(z[h]))
  grad@alpha_2h <- sapply(grad@alpha_2h, function(h) grad@alpha*weights@w_h[h]*f_d(z[h]))

  # Part 4
  grad@q_mh <- outer(1:M, 1:H , FUN=function(m,h) grad@alpha_2h[h]*f(r[m]) )

  # Part 5
  grad@alpha_1m <- sapply(grad@alpha_1m, function(m) grad@alpha*sum(weights@w_h*f_d(z)*weights@q_mh[m, ]*f_d(r[m])) )

  # Part 6
  grad@w_im <- outer(1:I, 1:M , FUN=function(i,m) grad@alpha_1m[m]*x )

  return(grad)
  
}# end function