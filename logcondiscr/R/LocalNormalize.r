LocalNormalize <- function(psi, dX){

# normalize psi s.t. the Lagrange term of L is 1

m <- length(psi)
last <- 0
if (m > 2){last <- sum(exp(psi[2:(m - 1)]))}
lagrangeTerm <- sum(J00(psi[1:(m - 1)], psi[2:m], dX)) - last
psi_new <- psi - log(lagrangeTerm)

return(psi_new)
}
