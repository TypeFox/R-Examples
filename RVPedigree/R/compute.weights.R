##' Compute/create matrix with genotype weights
##'
##' If the \code{weights} parameter (a numeric vector) is given, this
##'     function will return a matrix with that vector as the
##'     diagonal. If the \code{weights} parameter is not given (equal
##'     to \code{NULL}) the genotype weights will be calculated using
##'     the beta distribution (see the explanation of the
##'     \code{weights} parameter).
##' @title Compute/create matrix with genotype weights
##' @inheritParams RVPedigree
##' @param G matrix of genotypes (genetic variants as columns,
##'     individuals as rows) as output by the \code{\link{Get.G}}
##'     function.
##' @return Matrix with genotype weights.
##' @author Karim Oualkacha, Lennart C. Karssen
##' @keywords internal
compute.weights <- function(G, weights){
    if (is.null(weights)) {
        ## Calculate the minor allele frequency.
        freq.MAF <- apply(G, 2.0, mean)

        ## Set weights of the variants
        w <- (dbeta(freq.MAF, 1, 25))^2
        weights <- diag(w)
    } else {
        weights <- diag(weights)
    }

    return(weights)
}
