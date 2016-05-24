# Functions for MD ciriterion
# requires solve_LASP from package clue


# subfunction that does the minimization over a permutation matrix
pMatrix.min <- function(A){
   cost <- t(apply(A^2, 1, sum) - 2 * A + 1)
   vec <- c(solve_LSAP(cost))
   list(A=A[vec,], pvec=vec)
}

#pMatrix.min <- function(A, B) {
#    n <- nrow(A)
#    cost <- matrix(NA, n, n)
#    for (i in 1:n) {
#    for (j in 1:n) {
#        cost[j, i] <- (sum((B[j, ] - A[i, ])^2)) # correct Frobenius norm
#    } }
#    vec <- c(solve_LSAP(cost))
#    list(A=A[vec,], pvec=vec)
#    }

# main function
# input: square mixing matrix A
#        square unmixing matrix W.hat

MD <- function(W.hat,A)
    {
    G <- W.hat %*% A
    RowNorms <- sqrt(rowSums(G^2))
    G.0 <- sweep(G,1,RowNorms, "/")
    G.tilde <- G.0^2
    
    p <- nrow(A)
 #   B <- diag(p)
 #   Pmin <- pMatrix.min(G.tilde, B)
    Pmin <- pMatrix.min(G.tilde)
    G.tilde.p <- Pmin$A
    
    md <- sqrt(p - sum(diag(G.tilde.p)))/sqrt(p-1)
    md
    }
