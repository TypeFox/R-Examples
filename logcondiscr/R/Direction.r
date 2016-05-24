Direction <- function(W, psi, dX){

# Output:
# stepDirection : vector that points in direction of max L(psi) via Newton step

mat <- as.matrix(HesseL(psi, dX))
stepDirection <- - solve(mat) %*% GradientL(W, psi, dX)
return(stepDirection)
}
