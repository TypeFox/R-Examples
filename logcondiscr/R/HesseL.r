HesseL <- function(psi, dX){

# Compute Hesse matrix of log-likelihood L at psi

m <- length(psi)

# main diagonal
diag <- rep(0, m)
diag[1] <- -J20(psi[1], psi[2], dX[1])
diag[m] <- -J20(psi[m], psi[m - 1], dX[length(dX)]) 
if (m > 2){
    tmp1 <- J20(psi[2:(m - 1)], psi[1:(m - 2)], dX[1:(length(dX) - 1)])
    tmp2 <- J20(psi[2:(m - 1)], psi[3:m], dX[2:length(dX)])
    diag[2:(m - 1)] <- - tmp1 - tmp2 - exp(psi[2:(m - 1)])
    }

# off diagonals
tmpDiag <- -J11(psi[1:(m - 1)], psi[2:m], dX[1:length(dX)])
obDiag <- c(0, tmpDiag)
unDiag <- c(tmpDiag, 0)

# collect results in matrix
hesseL <- diag(m) * 0
hesseL <- Matrix::bandSparse(n = m, k = -1:1, diagonals = list(unDiag, diag, obDiag[-1]))
return(hesseL)
}
