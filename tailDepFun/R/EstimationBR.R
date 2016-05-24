#' @include helpFunctions.R
#' @include Other.R
NULL

#' Asymptotic variance matrix for the Brown-Resnick process.
#'
#' Computes the asymptotic variance matrix for the Brown-Resnick process, estimated using the pairwise M-estimator or the weighted least squares estimator.
#'
#' @param locations A \eqn{d} x 2 matrix containing the Cartesian coordinates of \eqn{d} points in the plane.
#' @param indices A \eqn{q} x \eqn{d} matrix containing exactly 2 ones per row, representing a pair of points from the matrix \code{locations}, and zeroes elsewhere.
#' @param par The parameters of the Brown-Resnick process. Either \eqn{(\alpha,\rho)} for an isotropic process or \eqn{(\alpha,\rho,\beta,c)} for an anisotropic process.
#' @param method Choose between "Mestimator" and "WLS".
#' @param Tol For "Mestimator" only. The tolerance in the numerical integration procedure. Defaults to 1e-05.
#' @return A \code{q} by \code{q} matrix.
#' @references Einmahl, J.H.J., Kiriliouk, A., Krajina, A., and Segers, J. (2016). An Mestimator of spatial tail dependence. Journal of the Royal Statistical Society: Series B (Statistical Methodology), 78(1), 275-298.
#' @references Einmahl, J.H.J., Kiriliouk, A., and Segers, J. (2016). A continuous updating weighted least squares estimator of tail dependence in high dimensions. See http://arxiv.org/abs/1601.04826.
#' @seealso \code{\link{selectGrid}}
#' @details The parameters of a The matrix \code{indices} can be either user-defined or returned from the function \code{selectGrid} with \code{cst = c(0,1)}. Calculation might be rather slow for \code{method = "Mestimator"}.
#' @export
#' @examples
#' locations <- cbind(rep(1:2, 3), rep(1:3, each = 2))
#' indices <- selectGrid(cst = c(0,1), d = 6, locations = locations, maxDistance = 1)
#' AsymVarBR(locations, indices, par = c(1.5,3), method = "WLS")
AsymVarBR<-function(locations, indices, par, method, Tol = 1e-05){
    if(length(par) == 2){
        if(any(par <= 0) || par[1] >= 2){warning("Invalid parameter values")}
    } else{
        if(any(par <= 0) || par[1] >= 2 || par[3] > pi/2){warning("Invalid parameter values")}
    }
    if(any(rowSums(indices) != 2)){
        warning("Only pairs are accepted: adjust indices accordingly")
     } else{
        coord <- coordinates(locations, indices)
        if(method == "Mestimator"){
            return(MestimatorAsymVarBR(pairs=coord,pars=par,Tol=Tol))
        } else if(method == "WLS"){
            return(WLSAsymVarBR(J=coord,pars=par))
        } else{
            warning("invalid method")
        }
    }
}

MestimatorBR <- function(x, locations, indices, k, isotropic, Tol, startingValue, Omega,iterate,covMat) {
    ranks <- apply(x, 2, rank)
    tuples<-lapply(seq(1:nrow(indices)), function(I) ranks[,indices[I,] == TRUE])
    totL<-unlist(lapply(tuples, function(i) stdfEmpInt(i,k=k)))
    pairs <- coordinates(locations, indices)
    if(isotropic){
        theta <- thetaPilot <- stats::optim(par = startingValue, fn = MestimatorMinimizeBR, method="BFGS", totlist = totL,
                                     pairs = pairs,w = Omega,control=list(maxit=10000,reltol=1e-12))
        if (iterate) {
            Omega <- solve(MestimatorAsymVarBR(pairs, thetaPilot$par, Tol))
            theta <- stats::optim(par = thetaPilot$par,fn = MestimatorMinimizeBR,method="BFGS",totlist = totL,
                           pairs = pairs,w = Omega,control=list(maxit=10000,reltol=1e-12))
        }
        covMatrix<-NULL
        if(covMat){
            psidot <- MestimatorpsidotBRIso(theta$par, pairs)
            temp <- solve(t(psidot) %*% Omega %*% psidot) %*% t(psidot)
            GammaMat <- MestimatorAsymVarBR(pairs, theta$par, Tol)
            covMatrix <- (temp %*% Omega %*% GammaMat %*% Omega %*% t(temp)) / k
        }
        return(list(theta = theta$par,
                    theta_pilot = thetaPilot$par,
                    covMatrix = covMatrix,
                    value = theta$value))

    } else{
        tau<-parsBRtoTau(startingValue[2:4])
        alpha<-startingValue[1]
        theta <- thetaPilot <- stats::optim(par = c(alpha,tau), fn = MestimatorMinimizeBR, method="BFGS", totlist = totL,
                                 pairs = pairs,w = Omega,control=list(maxit=10000,reltol=1e-12))
        if (iterate) {
          Omega <- solve(MestimatorAsymVarBR(pairs, thetaPilot, Tol))
          theta <- stats::optim(par = thetaPilot$par,fn = MestimatorMinimizeBR,method="BFGS",totlist = totL,
                       pairs = pairs,w = Omega,control=list(maxit=10000,reltol=1e-12))
        }
        covMatrix<-NULL
        if(covMat){
            psidot <- MestimatorpsidotBR(theta$par, pairs)
            temp <- solve(t(psidot) %*% Omega %*% psidot) %*% t(psidot)
            GammaMat <- MestimatorAsymVarBR(pairs, theta$par, Tol)
            covMatrixTau <- (temp %*% Omega %*% GammaMat %*% Omega %*% t(temp)) / k
            tauDeriv <- tauDerivMatrix(parsTauToBR(theta$par[2:4]))
            covMatrix <- t(tauDeriv) %*% covMatrixTau %*% tauDeriv
        }
        return(list(theta = c(theta$par[1],parsTauToBR(theta$par[2:4])),
                theta_pilot = c(thetaPilot$par[1],parsTauToBR(thetaPilot$par[2:4])),
                covMatrix = covMatrix,
                value = theta$value))
    }
}

WLSestimatorBR <- function(x, locations, indices, k, biascorr, k1, tau, isotropic, startingValue, Omega, iterate, covMat){
    ranks <- apply(x, 2, rank)
    tuples<-lapply(seq(1:nrow(indices)), function(I) ranks[,indices[I,] == TRUE])
    if(biascorr){
        totL<-unlist(lapply(tuples, function(i) stdfEmpCorr(i, k=k, tau = tau, k1 = k1)))
    } else{
        totL<-unlist(lapply(tuples, function(i) stdfEmp(i,k = k)))
    }
    J <- coordinates(locations, indices)
    if(isotropic){
        theta <- thetaPilot <- stats::optim(startingValue, fn = WLSminimizeBR,method = "L-BFGS-B", lower = c(0.01,0.01), upper= c(1.99,10),
                                     J = J, totlist = totL, w = Omega)
        if (iterate) {
            theta <- stats::optim(thetaPilot$par, fn = WLSminimizeBRcu, method = "L-BFGS-B", lower = c(0.01,0.01), upper= c(1.99,10),
                            J = J, totlist = totL)
        }
        covMatrix<-NULL
        if(covMat){
            phidot <- phidotBRECIso(J,theta$par)
            GammaMat <- WLSAsymVarBR(J,theta$par)
            if(iterate){
                covMatrix <- solve(t(phidot) %*% solve(GammaMat) %*% phidot) / k
            } else{
                temp <- solve(t(phidot) %*% Omega %*% phidot) %*% t(phidot)
                covMatrix <- (temp %*% Omega %*% GammaMat %*% Omega %*% t(temp)) / k
            }
        }

        return(list(theta = theta$par,
                    theta_pilot = thetaPilot$par,
                    covMatrix = covMatrix,
                    value = theta$value))
    } else{
        tau<-parsBRtoTau(startingValue[2:4])
        alpha<-startingValue[1]
        theta <- thetaPilot <- stats::optim(c(alpha,tau), fn = WLSminimizeBR,method = "BFGS",
                                     J = J, totlist = totL, w = Omega)
        if (iterate) {
            theta <- stats::optim(thetaPilot$par, fn = WLSminimizeBRcu, method = "BFGS",
                           J = J, totlist = totL)
        }
        covMatrix<-NULL
        if(covMat){
            phidot <- phidotBREC(J,theta$par)
            GammaMat <- WLSAsymVarBR(J,theta$par)
            if(iterate){
                covMatrixTau <- solve(t(phidot) %*% solve(GammaMat) %*% phidot) / k
            } else{
                temp <- solve(t(phidot) %*% Omega %*% phidot) %*% t(phidot)
                covMatrixTau <- (temp %*% Omega %*% GammaMat %*% Omega %*% t(temp)) / k
            }
            tauDeriv <- tauDerivMatrix(parsTauToBR(theta$par[2:4]))
            covMatrix <- t(tauDeriv) %*% covMatrixTau %*% tauDeriv
        }

        return(list(theta = c(theta$par[1],parsTauToBR(theta$par[2:4])),
                    theta_pilot = c(thetaPilot$par[1],parsTauToBR(thetaPilot$par[2:4])),
                    covMatrix = covMatrix,
                    value = theta$value))
    }
}


#' Estimation of the parameters of the Brown-Resnick process
#'
#' Estimation the parameters of the Brown-Resnick process, using either the pairwise M-estimator or weighted least squares (WLS).
#'
#' @param x An \eqn{n} x \eqn{d} data matrix.
#' @param locations A \eqn{d} x 2 matrix containing the Cartesian coordinates of \eqn{d} points in the plane.
#' @param indices A \eqn{q} x \eqn{d} matrix containing exactly 2 ones per row, representing a pair of points from the matrix \code{locations}, and zeroes elsewhere.
#' @param k An integer between 1 and \eqn{n - 1}; the threshold parameter in the definition of the empirical stable tail dependence function.
#' @param method Choose between \code{Mestimator} and \code{WLS}.
#' @param isotropic A Boolean variable. If \code{FALSE} (the default), then an anisotropic process is estimated.
#' @param Tol For \code{method = "Mestimator"} only. The tolerance parameter used in the numerical integration procedure. Defaults to 1e-05.
#' @param biascorr For \code{method = "WLS"} only. If \code{TRUE}, then the bias-corrected estimator of the stable tail dependence function is used. Defaults to \code{FALSE}.
#' @param k1 For \code{biascorr = TRUE} only. The value of \eqn{k_1} in the definition of the bias-corrected estimator of the stable tail dependence function.
#' @param tau For \code{biascorr = TRUE} only. The parameter of the power kernel.
#' @param startingValue Initial value of the parameters in the minimization routine. Defaults to \eqn{c(1,1.5)} if \code{isotropic = TRUE} and \eqn{c(1, 1.5, 0.75, 0.75)} if \code{isotropic = FALSE}.
#' @param Omega A \eqn{q} x \eqn{q} matrix specifying the metric with which the distance between the parametric and nonparametric estimates will be computed. The default is the identity matrix, i.e., the Euclidean metric.
#' @param iterate A Boolean variable. If \code{TRUE}, then for \code{method = "Mestimator"} the estimator is calculated twice, first with \code{Omega} specified by the user, and then a second time with the optimal \code{Omega} calculated at the initial estimate. If \code{method = "WLS"}, then continuous updating is used.
#' @param covMat A Boolean variable. If \code{TRUE} (the default), the covariance matrix is calculated. Standard errors are obtained by taking the square root of the diagonal elements.
#' @references Einmahl, J.H.J., Kiriliouk, A., and Segers, J. (2016). A continuous updating weighted least squares estimator of tail dependence in high dimensions. See http://arxiv.org/abs/1601.04826.
#' @references Einmahl, J.H.J., Kiriliouk, A., Krajina, A., and Segers, J. (2016). An Mestimator of spatial tail dependence. Journal of the Royal Statistical Society: Series B (Statistical Methodology), 78(1), 275-298.
#' @return A list with the following components:
#' \tabular{ll}{
#' \code{theta} \tab The estimator using the optimal weight matrix. \cr
#' \code{theta_pilot} \tab The estimator without the optimal weight matrix.\cr
#' \code{covMatrix} \tab The estimated covariance matrix for the estimator. \cr
#' \code{value} \tab The value of the minimized function at \code{theta}. \cr
#' }
#' @seealso \code{\link{selectGrid}}
#' @details The parameters of the Brown-Resnick process are either \eqn{(\alpha,\rho)} for an isotropic process or \eqn{(\alpha,\rho,\beta,c)} for an anisotropic process. The matrix \code{indices} can be either user-defined or returned from the function \code{selectGrid} with \code{cst = c(0,1)}. Estimation might be rather slow when \code{iterate = TRUE} or even when \code{covMat = TRUE}.
#' @export
#' @examples
#' ## define the locations of 9 stations
#' ## locations <- cbind(rep(c(1:3), each = 3), rep(1:3, 3))
#' ## select the pairs of locations
#' ## indices <- selectGrid(cst = c(0,1), d = 9, locations = locations, maxDistance = 1.5)
#' ## The Brown-Resnick process
#' ## set.seed(1)
#' ## x <- SpatialExtremes::rmaxstab(n = 1000, coord = locations, cov.mod = "brown",
#' ##                               range = 3, smooth = 1)
#' ## Calculate the estimtors.
#' ## EstimationBR(x, locations, indices, 100, method = "Mestimator", isotropic = TRUE,
#' ##             covMat = FALSE)$theta
#' ## EstimationBR(x, locations, indices, 100, method = "WLS", isotropic = TRUE,
#' ## covMat = FALSE)$theta

EstimationBR <- function(x, locations, indices, k, method, isotropic = FALSE, biascorr = FALSE,
                         Tol = 1e-05, k1 = (nrow(x) - 10), tau = 5, startingValue = NULL,
                       Omega = diag(nrow(indices)), iterate = FALSE, covMat=TRUE) {
    if(!is.null(startingValue)){
        if(any(startingValue < 1e-05) || startingValue[1] >= 2){
            warning("Incorrect startingValue: default is used")
            startingValue <- NULL
        }
        if(isotropic == FALSE){if(startingValue[3] < pi/2){
            warning("Incorrect startingValue: default is used")
            startingValue <- NULL
        }}
    } else{
        if(isotropic){
            startingValue <- c(1,1.5)
        } else{
            startingValue<-c(1,1.5,0.75,0.75)
        }
    }
    if(any(rowSums(indices) > 2) || (unique(c(indices)) != c(0,1) && unique(c(indices)) != c(1,0))){
        warning("The matrix indices should only contain rows with (d-2) zeroes and 2 ones")
    } else if(method == "Mestimator" && biascorr == TRUE){
        warning("biascorr can only be used for method = WLS")
    } else{
        if(method == "Mestimator"){
            return(MestimatorBR(x, locations, indices,k, isotropic = isotropic, Tol=Tol,
                    startingValue = startingValue,Omega=Omega,iterate=iterate,covMat=covMat))
        } else if(method == "WLS"){
            return(WLSestimatorBR(x, locations, indices, k, isotropic = isotropic,biascorr = biascorr,
                    k1 = k1, tau = tau, startingValue = startingValue, Omega = Omega,iterate=iterate,covMat=covMat))
        } else{
            warning("invalid method")
        }
    }
}
