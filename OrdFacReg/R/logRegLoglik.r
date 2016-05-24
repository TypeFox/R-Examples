logRegLoglik <- function(beta, Y, Z){

n <- length(Y)
pro <- Z %*% beta
L <- - sum(- pro * Y + log(1 + exp(pro)))

return(list("L" = L))
}
