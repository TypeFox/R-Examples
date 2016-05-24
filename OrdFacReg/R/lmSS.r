lmSS <- function(beta, D, Y){

L <- - sum((D - Y %*% beta) ^ 2)
dL <- 2 * t(Y)%*%((D - Y %*% beta))

return(list("L" = L, "dL" = dL))
}







