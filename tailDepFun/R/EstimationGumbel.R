#' @include helpFunctions.R
#' @include Other.R
NULL

#' Asymptotic variance matrix for the Gumbel model.
#'
#' Computes the asymptotic variance matrix for the Gumbel model, estimated using the pairwise M-estimator or the weighted least squares estimator.
#'
#' @param indices A \eqn{q} x \eqn{d} matrix containing at least 2 non-zero elements per row, representing the values in which we will evaluate the stable tail dependence function. For \code{method = Mestimator}, this matrix should contain exactly two ones per row.
#' @param par The parameter of the Gumbel model.
#' @param method Choose between "Mestimator" and "WLS".
#' @return A \code{q} by \code{q} matrix.
#' @seealso \code{\link{selectGrid}}
#' @references Einmahl, J.H.J., Kiriliouk, A., Krajina, A., and Segers, J. (2016). An Mestimator of spatial tail dependence. Journal of the Royal Statistical Society: Series B (Statistical Methodology), 78(1), 275-298.
#' @references Einmahl, J.H.J., Kiriliouk, A., and Segers, J. (2016). A continuous updating weighted least squares estimator of tail dependence in high dimensions. See http://arxiv.org/abs/1601.04826.
#' @details The matrix \code{indices} can be either user defines or returned by \code{selectGrid}. For \code{method = "Mestimator"}, only a grid with exactly two ones per row is accepted, representing the pairs to be used.
#' @export
#' @examples
#' indices <- selectGrid(c(0,1), d = 3, nonzero = c(2,3))
#' AsymVarGumbel(indices, par = 0.7, method = "WLS")
AsymVarGumbel<-function(indices, par, method){
    if(method == "Mestimator" && any(rowSums(indices) != 2)){
        warning("indices must contain exactly two ones per row for method = Mestimator")
    } else if(par <= 0 || par > 1){
        warning("par must be between zero and one")
    } else{
        if(method == "Mestimator"){
            return(AsymVarMestimatorGumbel(indices = indices, theta = par))
        } else if(method == "WLS"){
            return(AsymVarWLSGumbel(theta = par, indices = indices))
        } else{
            warning("invalid method")
        }
    }
}

WLSestimatorGumbel <- function(x, indices, k, biascorr, k1, tau, covMat) {
    ranks <- apply(x, 2, rank)
    if(biascorr){
        totL<-apply(indices, 1, function(j) stdfEmpCorr(ranks, k = k,tau = tau, k1 = k1, cst = j))
    } else{
        totL<-apply(indices, 1, function(j) stdfEmp(ranks, k = k, cst = j))
    }
    theta <- stats::optimize(f = WLSminimizeGumbel,interval = c(0.01,0.999), totlist = totL,
                      indices = indices)
    covMatrix<-NULL
    if(covMat){
        phidot <- psiWLSGumbel(theta$minimum, indices)
        GammaMat <- AsymVarWLSGumbel(indices = indices, theta = theta$minimum)
        temp <- solve(t(phidot) %*% phidot) %*% t(phidot)
        covMatrix <- (temp  %*% GammaMat %*%  t(temp)) / k
    }
    return(list(theta = theta$minimum, covMatrix = covMatrix, value = theta$objective))
}


MestimatorGumbel <- function(x, indices, k, covMat){
    q <- nrow(indices)
    ranks <- apply(x, 2, rank)
    tuples <- lapply(seq(1:q), function(I) ranks[,indices[I,] == TRUE])
    totL <- unlist(lapply(tuples, function(i) stdfEmpInt(i,k=k)))
    theta <- stats::optimize(f = MestimatorMinimizeGumbel, interval = c(0.01,0.999), q = q, totlist = totL)
    covMatrix<-NULL
    if(covMat){
        phidot <- psiMestimatorGumbel(theta$minimum, q)
        GammaMat <- AsymVarMestimatorGumbel(indices = indices, theta = theta$minimum)
        temp <- solve(t(phidot) %*% phidot) %*% t(phidot)
        covMatrix <- (temp  %*% GammaMat %*%  t(temp)) / k
    }
    return(list(theta = theta$minimum, covMatrix = covMatrix, value = theta$objective))
}


#' Estimation of the parameter of the Gumbel model
#'
#' Estimation the parameter of the Gumbel model, using either the pairwise M-estimator or weighted least squares (WLS).
#'
#' @param x An \eqn{n} x \eqn{d} data matrix.
#' @param indices A \eqn{q} x \eqn{d} matrix containing at least 2 non-zero elements per row, representing the values in which we will evaluate the stable tail dependence function. For \code{method = Mestimator}, this matrix should contain exactly two ones per row.
#' @param k An integer between 1 and \eqn{n - 1}; the threshold parameter in the definition of the empirical stable tail dependence function.
#' @param method Choose between \code{Mestimator} and \code{WLS}.
#' @param biascorr For \code{method = "WLS"} only. If \code{TRUE}, then the bias-corrected estimator of the stable tail dependence function is used. Defaults to \code{FALSE}.
#' @param k1 For \code{biascorr = TRUE} only. The value of \eqn{k_1} in the definition of the bias-corrected estimator of the stable tail dependence function.
#' @param tau For \code{biascorr = TRUE} only. The parameter of the power kernel.
#' @param covMat A Boolean variable. If \code{TRUE} (the default), the covariance matrix is calculated. Standard errors are obtained by taking the square root of the diagonal elements.
#' @details The matrix \code{indices} can be either user defined or returned by \code{selectGrid}. For \code{method = "Mestimator"}, only a grid with exactly two ones per row is accepted, representing the pairs to be used.
#' @return For \code{WLS}, a list with the following components:
#' \tabular{ll}{
#' \code{theta} \tab The estimator with weight matrix identity. \cr
#' \code{covMatrix} \tab The estimated covariance matrix for the estimator. \cr
#' \code{value} \tab The value of the minimized function at \code{theta}. \cr
#' }
#' @seealso \code{\link{selectGrid}}
#' @references Einmahl, J.H.J., Kiriliouk, A., and Segers, J. (2016). A continuous updating weighted least squares estimator of tail dependence in high dimensions. See http://arxiv.org/abs/1601.04826.
#' @references Einmahl, J.H.J., Kiriliouk, A., Krajina, A., and Segers, J. (2016). An Mestimator of spatial tail dependence. Journal of the Royal Statistical Society: Series B (Statistical Methodology), 78(1), 275-298.
#' @export
#' @examples
#' ## Generate data with theta = 0.5
#' ## set.seed(1)
#' ## n <- 1000
#' ## cop <- copula::gumbelCopula(param = 2, dim = 3)
#' ## data <- copula::rCopula(n = n,copula = cop)
#' ## Transform data to unit Pareto margins
#' ## x <- apply(data, 2, function(i) n/(n + 0.5 - rank(i)))
#' ## Define indices in which we evaluate the estimator
#' ## indices <- selectGrid(c(0,1), d = 3)
#' ## EstimationGumbel(x, indices, k = 50, method = "WLS", biascorr = TRUE)
EstimationGumbel <- function(x, indices, k, method, biascorr = FALSE, k1 = (nrow(x) - 10),tau = 5,covMat=TRUE) {
    if(method == "Mestimator" && biascorr == TRUE){
      warning("biascorr can only be used for method = WLS")
    } else{
    if(method == "Mestimator"){
        if(any(rowSums(indices) != 2)){
            warning("For method = Mestimator, only pairs are accepted: adjust indices accordingly")
        } else{
            return(MestimatorGumbel(x, indices = indices, k = k, covMat = covMat))
        }
     } else if(method == "WLS"){
            return(WLSestimatorGumbel(x, indices = indices, k = k, biascorr = biascorr, k1 = k1, tau = tau,covMat=covMat))
        } else{
            warning("invalid method")
        }
    }
}
