rIG <-
function(epsilon, Tau) {
## generating random number from inverse-gaussian dist'n
     z <- rchisq(1, 1)
     ss <- sqrt(4 * epsilon * z / Tau + (epsilon * z)^2)
     z1 <- epsilon + (epsilon^2) * Tau * z / 2 - (epsilon * Tau / 2) * ss
     u1 <- runif(1, 0, 1)
     xxx <- z1
     if(u1 > (epsilon / (epsilon + z1)) ) {
           xxx <- (epsilon^2) / z1
     }
     return(as.numeric(xxx))      
}
