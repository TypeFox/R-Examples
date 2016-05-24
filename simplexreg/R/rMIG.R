rMIG <-
function(epsilon, Tau, mu) {
## generating random number from inverse-gaussian mixture dist'n
     x1 <- rIG(epsilon, Tau)
     x2 <- rchisq(1, 1)     
     x3 <- x2 * Tau * (epsilon^2)
     u2 <- runif(1, 0, 1)
     xx <- x1
     if(u2 < mu) { xx <- x1 + x3 }
     return(as.numeric(xx))
}
