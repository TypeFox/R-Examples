#' Compute B and S matrices in algorithm of Vardi (1996)
#'
#' Function to compute B and S matrices for moment equations of Vardi's method
#' (1996). It's not particularly efficient, but it works.
#'
#' @param A routing matrix (m x k)
#' @param Y matrix of link loads over time (m x n, one column per time)
#' @return list containing two entries for the B and S matrices, respectively
#' @references Y. Vardi. Network tomography: estimating source-destination
#' traffic intensities from link data. Journal of the American Statistical
#' Association, 91:365-377, 1996.
#' @export
#' @keywords models multivariate
#' @family vardi
vardi.compute.BS <- function(A,Y) {
    r <- dim(A)[1]
    K <- dim(Y)[2]
    first <- 1
    for (i in 1:r) {
        for (j in i:r) {
            bi <- A[i,]*A[j,]
            si <- 1/K * sum(Y[i,]*Y[j,]) - mean(Y[i,])*mean(Y[j,])
            # Delete the impossible moment equations under the Poisson model.
            if (si < 0 && sum(bi) > 0 || si > 0 && sum(bi) < 0 ||
                si != 0 && sum(bi) == 0 || si == 0 && sum(bi) != 0) {
                next
            }
            if (first) {
                B <- bi
                S <- si
                first <- 0
                next
            }
            B <- rbind(B,bi)
            S <- c(S,si)
        }
    }
    dimnames(B) <- NULL

    return(list(B=B,S=S))
}

#' Execute single iteration for algorithm of Vardi (1996)
#'
#' Function to compute B and S matrices for moment equations of Vardi's method
#' (1996). It's not particularly efficient, but it works.
#'
#' @param A routing matrix (m x k)
#' @param yBar numeric vector of mean link loads (length m)
#' @param lambda value of lambda from last iteration
#' @param B B matrix computed by \code{\link{vardi.compute.BS}}
#' @param S S matrix computed by \code{\link{vardi.compute.BS}}
#' @return numeric vector of length k with updated lambda 
#' @references Y. Vardi. Network tomography: estimating source-destination
#' traffic intensities from link data. Journal of the American Statistical
#' Association, 91:365-377, 1996.
#' @export
#' @keywords models multivariate
#' @family vardi
vardi.iteration <- function(A, yBar, lambda, B, S) {
    adot <- colSums(A)
    bdot <- colSums(B)
    Alambda <- as.vector(crossprod(t(A),lambda))    # A %*% lambda
    Blambda <- as.vector(crossprod(t(B),lambda))    # B %*% lambda
    Ay <- A*yBar
    BS <- B*S
    AyoAlambda <- Ay/Alambda
    BSoBlambda <- BS/Blambda
    lambda.hat.Ayl <- lambda/apply(A,2,sum)*apply(AyoAlambda,2,sum)
    lambda.hat.BSl <- lambda/apply(B,2,sum)*apply(BSoBlambda,2,sum)

    lambda.new <- adot/(adot+bdot)*lambda.hat.Ayl +
    bdot/(adot+bdot)*lambda.hat.BSl

    return(lambda.new)
}

#' Run algorithm of Vardi (1996) given B and S matrices
#'
#' Runs moment-matching algorithm of Vardi (1996) until convergence
#'
#' @param A routing matrix (m x k)
#' @param Y matrix of link loads over time (m x n, one column per time)
#' @param lambda numeric vector of starting values for OD flows (length k)
#' @param tol numeric tolerance for halting iterations
#' @return numeric vector of length k with estimated OD flows
#' @references Y. Vardi. Network tomography: estimating source-destination
#' traffic intensities from link data. Journal of the American Statistical
#' Association, 91:365-377, 1996.
#' @export
#' @keywords models multivariate
#' @family vardi
vardi.algorithm <- function(A, Y, lambda, tol=1e-3) {
    # Compute B and S matrices
    BS <- vardi.compute.BS(A, Y)

    # Compute row means of Y
    yBar <- rowMeans(Y)

    # Run iterations until convergence
    lambda.new <- vardi.iteration(A,yBar,lambda,BS$B,BS$S)
    while (max(abs(lambda.new - lambda)) > tol) {
        lambda <- lambda.new
        lambda.new <- vardi.iteration(A,yBar,lambda,BS$B,BS$S)
    }

    return(lambda.new)
}

