##' @title Community estimation in G-models via CORD
##'
##' @description Partition data points (variables) into clusters/communities. Reference: Bunea, F., Giraud, C., & Luo, X. (2015). Community estimation in \eqn{G}-models via CORD. arXiv preprint arXiv:1508.01939. \url{http://arxiv.org/abs/1508.01939}.
##'
##' 
##' @param X Input data matrix. It should be an n (samples) by p (variables) matrix  when \code{input} is set to the value "data" by default. It can also be a p by p symmetric matrix when \code{X} is a correlation matrix or a distance matrix if \code{input} is set accordingly.
##' @param tau Threshold to use at each iteration. A theoretical choice is about \eqn{2n^{-1/2}\log^{1/2} p}.
##' @param kendall Whether to compute Kendall's tau correlation matrix from \code{X}, when \code{input} is set to "data". If \code{FALSE}, Pearson's correlation will be computed, usually faster for large p.
##' @param input Type of input \code{X}. It should be set to "data" when \code{X} is an  n (samples) by p (variables) matrix. If \code{X} is a correlation matrix or a distance matrix, it should be set to "cor" or "dist" respectively.
##' @return \code{list} with one element: a vector of integers showing which cluster/community each point is assigned to.
##'
##' @importFrom stats cor
##' 
##' @examples
##' set.seed(100)
##' X <- 2*matrix(rnorm(200*2), 200, 10)+matrix(rnorm(200*10), 200, 10)
##' cord(X)
##'
##' @keywords multivariate cluster
##' 
##' @export 
cord <- function(X, tau=2*sqrt(log(ncol(X))/nrow(X)), kendall=T,  input=c("data", "cor", "dist")) {
    
    input <- match.arg(input,c("data", "cor", "dist"))

    if(input == "cor") {
        input.cordDist <- F
        input.cor <- T
    } else if (input == "cordDist") {
        input.cordDist <- T
        input.cor <- F
    } else {
        input.cordDist <- F
        input.cor <- F
    }
     n <- nrow(X)
    p <- ncol(X)

    mat <- NULL
    
    
    if (input.cordDist) {
        fd <- X 
    } else {
        if (input.cor) {
            mat <- X
        } else {
            if (is.null(mat)) {
                ## compute distance matrix
                if (kendall) {
                    ## mat <- sin(pi*cor.fk(X)/2)
                    mat <- sin(pi*cor(X, method= "kendall")/2)
                } else {
                    mat <- cor(X) 
                }
            }
        }
        fd <- cordDist(mat)
    }
    
    ## add large values to the diagonal
    fd <- fd+diag(1e+5, p) 
    
    Gvec <- greedyCord(fd, tau)
    
    return(list(cluster=as.numeric(Gvec)))
}


