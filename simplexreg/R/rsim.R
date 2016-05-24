rsim <-
function(n, mu, sig) {
## generating random number from simplex dist'n 
## by transformation from inverse-gaussian mixture dist'n
     Sigma <- sig^2
     epsilon <- mu / (1 - mu)
     Tau <- Sigma * ((1 - mu)^2)
     yy <- rep(0, n)
     for(i in 1:n) {
              x <- rMIG(epsilon, Tau, mu)
              yy[i] <- x / (1 + x)
     }
     return(as.vector(yy))
}
