#' Matrix Of Eigenvalues Of A List Of Symetric Matrices
#' 
#' @description Computation of the eigenvalues of each matrix Vj in the case of two random effects (random =c(1,2)), done via \code{eigen}
#' @param V list of matrices Vj
#' @return
#' \item{eigenvalues}{Matrix of 2 rows and as much columns as matrices V}
#' @references See Bidimensional random effect estimation in mixed stochastic differential model, C. Dion and V. Genon-Catalot,  \emph{Stochastic Inference for Stochastic Processes 2015, Springer Netherlands}, \bold{1--28}





eigenvaluesV <- function(V) {
    
    eig <- matrix(0, length(V), 2)
    for (j in seq_along(V)) {
        eig[j, ] <- c(eigen(V[[j]], symmetric = TRUE)$values)
    }
    return(eigenvalues = eig)
} 
