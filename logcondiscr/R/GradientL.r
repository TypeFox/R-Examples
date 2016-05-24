GradientL <- function(W, psi, dX){

# Computes gradient of the log-likelihood L at psi
m <- length(psi)
tmp <- rep(NA, m)
tmp[1] <- J10(psi[1], psi[2], dX[1])
tmp[m] <- J10(psi[length(psi)], psi[length(psi) - 1], dX[length(dX)])

if (m > 2){
    J <- J10(psi[2:(m - 1)], psi[3:m], dX[2:length(dX)]) + J10(psi[2:(m - 1)], psi[1:(m - 2)], dX[1:(length(dX) - 1)]) - exp(psi[2:(m - 1)])
    tmp[2:(m - 1)] <- J
    }

gradL <- W - tmp
return(gradL)
}
