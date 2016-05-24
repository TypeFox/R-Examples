#' genComplResid
#'
#' Generates a completed vector of residuals
#'
#' This function involves three steps. The first two are similar in spirit to
#' the two-stage procedure of Othus and Li (2010).

#' \enumerate{
#' \item The vector of covariate parameters and the monotone increasing function
#' of the transformation model with censored data (Cheng et al., 1995) are 
#' estimated under the working independence assumption following the algorithm 
#' of Chen et al. (2002) and used to compute raw residuals;
#' \item The polygenic heritability parameter is estimated which is a measure of
#' the dependence between the survival traits of correlated groups that cannot 
#' be attributed to the SNP set under investigation. This estimate is used to 
#' deduce the approximate covariance matrix of the raw residuals.
#' \item An imputation procedure is employed to replace the censored raw 
#' residuals by the mean of multiple imputed values generated from the posterior
#' distribution of the uncensored version with the restriction to be larger than
#' the original censored values, componentwise. The completed vector of 
#' residuals is then deduced and standardized. A scale parameter is used to 
#' reflect the fact that we are using multiple imputed values rather than real 
#' observations.
#' }
#' \bold{Warning:} Correlated groups identified by the vector \code{blkID} most
#' often corresponds to families or blocks of the block-diagonal kinship matrix
#' \bold{Phi}. Larger groups such as regions of residence can be considered, for
#' example to take into account population stratification or cryptic relatedness.
#' However, the number of censored individuals in each group cannot exceed 1000
#' as the test makes use of the distribution function of the multivariate normal 
#' distribution for which the maximum dimension is 1000 in the function 
#' \bold{pmvnorm} of the package \pkg{mvtnorm}.
#'
#' Simulation studies reported in Leclerc et al. (2015) suggest that the use
#' of \code{m = 50} imputations guarantees a reasonable power in practice.
#'
#' \bold{Warning:} No missing data is allowed for \code{U}, \code{Delta}, 
#' \code{Phi}, \code{blkID}, and \code{X}.
#'
#' @param U a nx1 vector containing the survival times. \code{U = min(C, T)} 
#' where \code{C} is the censoring time, and \code{T} the failure time
#' @param Delta a nx1 vector containing the censoring indicator
#' @param Phi a nxn kinship matrix
#' @param blkID a nx1 vector with entries identifying correlated groups of 
#' observations. The number of censored individuals in each group cannot exceed
#' 1000 (see \bold{Details})
#' @param m default=50. Number of imputations used to generate the completed 
#' vector of residuals
#' @param X a nxp \bold{matrix} of p covariates. Each row represents a different 
#' individual, and each column represents a different numeric covariate. If no 
#' covariates are present, \code{X} can be left as \code{NULL}
#' @return The function produces a list consisting of:
#' @return \item{compResid}{the completed vector of residuals}
#' @return \item{herit}{the estimate of the polygenic heritability parameter}
#' @return \item{covPar}{the estimate of the vector of covariate parameters (if 
#' applicable)}
#' @author Martin Leclerc <martin.leclerc.5@@ulaval.ca> and Lajmi Lakhal Chaieb
#'  <lakhal@@mat.ulaval.ca>
#' @references Chen K, Jin Z, Ying Z. 2002. Semiparametric analysis of 
#' transformation models with censored data. Biometrika 89:659-668.
#' 
#' Cheng SC, Wei LJ, Ying Z. 1995. Analysis of transformation models with 
#' censored data. Biometrika 82:835-845.
#' 
#' Leclerc M, The Consortium of Investigators of Modifiers of BRCA1/2, Simard J,
#' Lakhal-Chaieb L. 2015. SNP set association testing for survival outcomes
#' in the presence of intrafamilial correlation. Genetic Epidemiology 
#' 39:406-414.
#' 
#' Othus M, Li Y. 2010. A gaussian copula model for multivariate survival data.
#' Stat Biosci 2:154-179.
#' @examples
#' data(simGyriq)
#' for (i in seq_along(simGyriq)) assign(names(simGyriq)[i], simGyriq[[i]])
#'
#' cr <- genComplResid(U, Delta, Phi, blkID, m=50, X)
#' @importFrom mvtnorm dmvnorm pmvnorm rmvnorm 
#' @export
genComplResid <- function(U, Delta, Phi, blkID, m=50, X=NULL) {

    if (class(Phi) != "matrix") stop("Phi is not a matrix")
    if (nrow(Phi) != ncol(Phi)) 
        stop("Phi is not a square matrix")
    
    if (sum(is.na(U)) != 0) stop("U cannot have any missing values")
    if (sum(is.na(Delta)) != 0) stop("Delta cannot have any missing values")
    if (sum(is.na(Phi)) != 0) stop("Phi cannot have any missing values")
    if (sum(is.na(blkID)) != 0) stop("blkID cannot have any missing values")
    
    if (length(U) != length(Delta)) 
        stop("Dimensions of U and Delta do not match")
    if (length(U) != nrow(Phi)) 
        stop("Dimensions of U and Phi do not match")
    if (length(U) != length(blkID)) 
        stop("Dimensions of U and blkID do not match")

    blkCens <- aggregate(x = 1-Delta, by = list(blkID), FUN = sum)$x
    if (max(blkCens) > 1000)
        stop("The number of censored individuals in each correlated group of
             observations identified by the vector 'blkID' cannot exceed 1000.")
        
    if (is.null(X) == TRUE) {
    	covPar <- NULL
        covTerm <- rep(0, length(U))
    } else {
        if (class(X) != "matrix") stop("X is not a matrix")
        if (length(U) != nrow(X)) stop("Dimensions of U and X do not match")
        
        if (sum(apply(X, 2, chkVariation)) != ncol(X)) 
            stop("One of the covariates in X has no variation")
        
        fit <- survival::coxph(survival::Surv(U, Delta) ~ X)
        covPar <- fit$coef
        covTerm <- as.vector(X %*% covPar)
    }

    rawResid <- estimRawResid(U=U, Delta=Delta, covTerm=covTerm)
    
    herit <- estimHeritPar(rawResid=rawResid, Delta=Delta, Phi=Phi, blkID=blkID)
    Gam <- herit * Phi + (1 - herit) * diag(rep(1, dim(Phi)[1]))
    
    if (sum(Delta) == length(Delta)) {
        GamEig <- eigen(Gam)
        GamEigVal <- GamEig$values
        GamSqrt <- GamEig$vectors %*% diag(sqrt(GamEigVal)) %*% 
            solve(GamEig$vectors)
        compResid <- as.vector(solve(GamSqrt) %*% rawResid)
    } else {
        compResid <- imputeResid(rawResid=rawResid, Delta=Delta, Gam=Gam, m=m, 
                                 blkID=blkID)
    }
    return(list(compResid=compResid, herit=herit, covPar=covPar))
}
