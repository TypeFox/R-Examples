ResExpG <-
function(r) { zero <- 1e-6
num <- dnorm(r); den <- 1-pnorm(r)
val <- num/den; ind <- (den < zero) ; val[ind] <- r[ind] 
val}

