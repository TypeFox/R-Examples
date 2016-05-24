logRegDeriv <- function(beta, Y, Z){

n <- dim(Z)[1]
p <- dim(Z)[2]

pro <- Z %*% beta
p <- exp(pro) / (1 + exp(pro))

G <- t(Z) %*% (p - Y)

return(list("dL" = -G))
}
