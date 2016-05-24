.weightStrucFuncMatrix <- function(coeff) {
  function(g) {
    stopifnot(.validateGraph(g))
    dist <- distanceMatrix(g)
    pis <- infoTheoreticGCM(g,  coeff = coeff, infofunct = "sphere")$pis
    len <- length(pis)

    # calculate diff[i, j] = abs(pis[i] - pis[j])
    diff <- abs(matrix(pis, len, len, TRUE) - matrix(pis, len, len, FALSE))

    polyMat <- 1 - diff / 2^dist
    diag(polyMat) <- 1
    polyMat
  }
}

weightStrucFuncMatrix_lin <- .weightStrucFuncMatrix("lin")
weightStrucFuncMatrix_exp <- .weightStrucFuncMatrix("exp")
