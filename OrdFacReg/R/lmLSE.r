lmLSE <- function(D, Y){

beta <- (MASS::ginv(t(Y) %*% Y)) %*% t(Y) %*% D
L <- lmSS(beta, D, Y)$L
return(list("L" = L, "beta" = beta))
}







