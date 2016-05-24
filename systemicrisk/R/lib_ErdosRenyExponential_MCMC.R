#### from devtools::use_rcpp()
#' @useDynLib systemicrisk
#' @importFrom Rcpp sourceCpp
#' @importFrom lpSolve lp
#' @import stats utils

###

#' @title Sample from the ERE model with given row and column sums
#'
#' @description
#' Samples from the Erdos Reny model with Exponential weights and
#' known marginals.  Runs a Gibbs sampler to do this. A starting
#' liabilities is generated via \code{\link{getfeasibleMatr}} before
#' \code{\link{steps_ERE}} is called.
#'
#'
#' @param l vector of interbank libabilities
#' @param a vector of interbank assets
#' @param p probability of existence of a link (either a numerical
#' value or a matrix). A single numerical value is converted into a
#' matrix with 0s on the diagonal.
#' @param lambda instensity parameters - either a numerical value or a
#' matrix with positive entries)
#' @inheritParams steps_ERE
#' @return List of simulation results
#'
#'@examples
#' l <- c(1,2.5,3)
#' a <- c(0.7,2.7,3.1)
#' L <- sample_ERE(l,a,p=0.5,lambda=0.25,nsamples=5,thin=20,burnin=10)
#' L
#'
#'
#' @export
sample_ERE <- function(l,
                                        a,
                                        p,
                                        lambda,
                                        nsamples=1e4,
                                        thin=1e3,
                                        burnin=1e4
                                        ){
    n <- length(l)
    if (!is.matrix(p)){
        p <- matrix(p,nrow=n,ncol=n);
        diag(p) <- 0
    }
    if (!is.matrix(lambda)){;
        lambda <- matrix(lambda,nrow=n,ncol=n);
        diag(lambda) <- 0
    }
    L <- findFeasibleMatrix(l,a,p)

    steps_ERE(L=L,p=p,lambda=lambda,nsamples=nsamples,thin=thin,burnin=burnin)
}

#' @title Perform Steps of the Gibbs Sampler of the ERE model
#'
#' @description
#' Runs a Gibbs sampler in the Erdos Reny model with Exponential weights (ERE model)
#' and fixed marginals. The algorithm starts from a given matrix.
#'
#' @param L Starting matrix for the Gibbs sampler. Implicitly defines the fixed marginals.
#' @param p A matrix with entries in [0,1]
#' @param lambda A matrix with nonnegative entries
#' @param nsamples Number of samples to return.
#' @param thin Frequency at which samples should be generated (default=1, every step)
#' @param burnin Number of initial steps to discard.
#' @return List of simulation results
#'
#' @examples
#' L <- matrix(rexp(4*4),nrow=4,ncol=4); diag(L)=0;
#' p <- matrix(0.5,nrow=4,ncol=4); diag(p) <-0;
#' lambda <- matrix(1,nrow=4,ncol=4); diag(lambda)<-0;
#'
#' L <- steps_ERE(L=L,p=p,lambda=lambda,nsamples=5,thin=50,burnin=20)
#' L
#'
#' @seealso \code{\link{sample_ERE}}
#'
#' @export
steps_ERE <- function(L,
                      p,
                      lambda,
                      nsamples=1e4,
                      thin=1e3,
                      burnin=1e4
                      ){
    L <- cloneMatrix(L)## to ensure that there are no strange side effects

    res <- list()
        GibbsSteps_kcycle(L=L,p=p,lambda=lambda,it=burnin)
        for (i in 1:nsamples){
            GibbsSteps_kcycle(L=L,p=p,lambda=lambda,it=thin)
            res[[i]] <- cloneMatrix(L)
        }
    res
}
