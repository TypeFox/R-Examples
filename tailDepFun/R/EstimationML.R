#' @include helpFunctions.R
#' @include Other.R
NULL

#' Asymptotic variance matrix for the max-linear model.
#'
#' Computes the asymptotic variance matrix for the max-linear model, estimated using the weighted least squares estimator.
#'
#' @param indices A \eqn{q} x \eqn{d} matrix containing at least 2 non-zero elements per row, representing the values in which we will evaluate the stable tail dependence function.
#' @param par The parameter vector.
#' @param Bmatrix A function that converts the parameter vector theta to a parameter matrix B. If \code{NULL}, then a simple 2-factor model is assumed.
#' @return A \code{q} by \code{q} matrix.
#' @seealso \code{\link{selectGrid}}
#' @references Einmahl, J.H.J., Kiriliouk, A., and Segers, J. (2016). A continuous updating weighted least squares estimator of tail dependence in high dimensions. See http://arxiv.org/abs/1601.04826.
#' @export
#' @examples
#' indices <- selectGrid(c(0,0.5,1), d = 3, nonzero = 3)
#' AsymVarMaxLinear(indices, par = c(0.1,0.55,0.75))
AsymVarMaxLinear<-function(indices, par, Bmatrix = NULL){
    if(is.null(Bmatrix)){
        Bpars <- cbind(par, 1 - par)
    } else{
        Bpars <- Bmatrix(par)
    }
    if(any(rowSums(Bpars) > 1) || any(c(Bpars) < 0)){
        warning("invalid parameter values")
    } else{
        return(AsymVarMLwls(Bpars=Bpars, indices = indices, d = ncol(indices)))
    }
}

WLSestimatorML <- function(x, indices, k, Bmatrix, phidot, biascorr, k1, tau, startingValue, Omega, iterate, covMat, GoFtest, dist) {
    ranks <- apply(x, 2, rank)
    if(biascorr){
        totL<-apply(indices, 1, function(j) stdfEmpCorr(ranks, k = k,tau = tau, k1 = k1, cst = j))
    } else{
        totL<-apply(indices, 1, function(j) stdfEmp(ranks, k = k, cst = j))
    }
    theta <- thetaPilot <- stats::optim(startingValue, fn = WLSminimizeML,totlist = totL, indices = indices,
                                 Bmatrix = Bmatrix, w = Omega, control=list(maxit=20000,reltol=1e-10))
    if(iterate){
        theta <- stats::optim(thetaPilot$par, fn = WLSminimizeMLcu,totlist = totL, indices = indices,
                       Bmatrix = Bmatrix, control=list(maxit=20000,reltol=1e-10))
    }
    covMatrix <- GammaMat <- GoFresult <- NULL
        if(covMat){
            phidot <- phidot(indices,theta$par)
            GammaMat <- AsymVarMLwls(Bpars = Bmatrix(theta$par), indices = indices, d = ncol(indices))
            if(iterate){
                while(rcond(GammaMat) < 1e-05){GammaMat <- GammaMat + (1e-05)*diag(length(totL))}
                covMatrix <- solve(t(phidot) %*% solve(GammaMat) %*% phidot) / k
            } else{
                temp <- solve(t(phidot) %*% Omega %*% phidot) %*% t(phidot)
                covMatrix <- (temp %*% Omega %*% GammaMat %*% Omega %*% t(temp)) / k
            }
            colnames(covMatrix) <- NULL
            rownames(covMatrix) <- NULL
        }
        if(GoFtest){
            prphi <- solve(t(phidot) %*% phidot)
            P <- phidot %*% prphi %*% t(phidot)
            if(is.null(GammaMat)){
                GammaMat <- AsymVarMLwls(Bpars = Bmatrix(theta$par), indices = indices, d = ncol(indices))
            }
            TotMat <- (diag(rep(1,nrow(indices))) - P)%*%GammaMat%*%t(diag(rep(1,nrow(indices))) - P)
            eig <- eigen(TotMat, symmetric = TRUE)
            s <- length(which(eig$values > dist))
            if(s == 1){
                A <- as.matrix(eig$vectors[,1:s]) %*% as.matrix(1/eig$values[1:s]) %*% t(as.matrix(eig$vectors[,1:s]))
            } else{
                A <- eig$vectors[,1:s] %*% diag(1/eig$values[1:s]) %*% t(eig$vectors[,1:s])
            }
            res <- WLSminimizeML(theta$par, totlist = totL, indices = indices,
                                 Bmatrix = Bmatrix, w = A)
            GoFresult = list(value = c(res)*k, s = s)
        }
        return(list(theta = theta$par,
                    theta_pilot = thetaPilot$par,
                    covMatrix = covMatrix,
                    value = theta$value,
                    GoFresult = GoFresult))
}


MestimatorML <- function(x, indices, k, Bmatrix, startingValue){
    ranks <- apply(x, 2, rank)
    mpairs <- lapply(1:nrow(indices), function(I) ranks[,indices[I,] == TRUE])
    totL<-unlist(lapply(mpairs, function(i) stdfEmpInt(i, k = k)))
    temp <- stats::optim(startingValue, fn = MestimatorMinimizeML, Bmatrix = Bmatrix,
            indices = indices, totlist = totL,control=list(reltol=1e-12))
    return(theta = temp$par)
}


#' Estimation of the parameters of the max-linear model
#'
#' Estimation the parameters of the max-linear model, using either the pairwise M-estimator or weighted least squares (WLS).
#'
#' @param x An \eqn{n} x \eqn{d} data matrix.
#' @param indices A \eqn{q} x \eqn{d} matrix containing at least 2 non-zero elements per row, representing the values in which we will evaluate the stable tail dependence function.
#' @param k An integer between 1 and \eqn{n - 1}; the threshold parameter in the definition of the empirical stable tail dependence function.
#' @param method Choose between \code{Mestimator} and \code{WLS}.
#' @param Bmatrix A function that converts the parameter vector theta to a parameter matrix B. If nothing is provided, then a simple 2-factor model is assumed.
#' @param Ldot For \code{method = "WLS"} only. A \eqn{q} x \eqn{p} matrix, where \eqn{p} is the number of parameters of the model. Represents the total derivative of the function L defined in Einmahl et al. (2016). If nothing is provided, then a simple 2-factor model is assumed.
#' @param biascorr For \code{method = "WLS"} only. If \code{TRUE}, then the bias-corrected estimator of the stable tail dependence function is used. Defaults to \code{FALSE}.
#' @param k1 For \code{biascorr = TRUE} only. The value of \eqn{k_1} in the definition of the bias-corrected estimator of the stable tail dependence function.
#' @param tau For \code{biascorr = TRUE} only. The parameter of the power kernel.
#' @param startingValue Initial value of the parameters in the minimization routine.
#' @param Omega A \eqn{q} x \eqn{q} matrix specifying the metric with which the distance between the parametric and nonparametric estimates will be computed. The default is the identity matrix, i.e., the Euclidean metric.
#' @param iterate A Boolean variable. For \code{method = "WLS"} only. If \code{TRUE}, then continuous updating is used. Defaults to \code{FALSE}.
#' @param covMat A Boolean variable. For \code{method = "WLS"} only. If \code{TRUE} (the default), the covariance matrix is calculated.
#' @param GoFtest A Boolean variable. For \code{method = "WLS"} only. If \code{TRUE}, then the goodness-of-fit test of Corollary 2.6 from Einmahl et al. (2016) is performed. Defaults to \code{FALSE}.
#' @param dist A positive scalar. If \code{GoFtest = TRUE}, only eigenvalues \eqn{\nu} larger than \code{dist} are used; see Corollary 2.6 (Einmahl et al., 2016). Defaults to 0.01.
#' @param EURO A Boolean variable. If \code{TRUE}, then the model from Einmahl et al. (2016, Section 4) is assumed, and the corresponding \code{Bmatrix} and \code{Ldot} are used.
#' @details The matrix \code{indices} can be either user defined or returned by \code{selectGrid}.
#' For \code{method = "Mestimator"}, only a grid with exactly two ones per row is accepted,
#' representing the pairs to be used. The functions \code{Bmatrix} and \code{Ldot} can be defined
#' such that they represent a max-linear model on a directed acyclic graph: see the vignette for this package for an example.
#' @return For \code{Mestimator}, the estimator \code{theta} is returned. For \code{WLS}, a list with the following components:
#' \tabular{ll}{
#' \code{theta} \tab The estimator with estimated optimal weight matrix. \cr
#' \code{theta_pilot} \tab The estimator without the optimal weight matrix.\cr
#' \code{covMatrix} \tab The estimated covariance matrix for the estimator. \cr
#' \code{value} \tab The value of the minimized function at \code{theta}. \cr
#' \code{GoFresult} \tab A list of length two, returning the value of the test statistic and \code{s}. \cr
#' }
#' @seealso \code{\link{selectGrid}}
#' @references Einmahl, J.H.J., Kiriliouk, A., and Segers, J. (2016). A continuous updating weighted least squares estimator of tail dependence in high dimensions. See http://arxiv.org/abs/1601.04826.
#' @export
#' @examples
#' ## Generate data
#' set.seed(1)
#' n <- 1000
#' fr <- matrix(-1/log(runif(2*n)), nrow = n, ncol = 2)
#' data <- cbind(pmax(0.3*fr[,1],0.7*fr[,2]),pmax(0.5*fr[,1],0.5*fr[,2]),pmax(0.9*fr[,1],0.1*fr[,2]))
#' ## Transform data to unit Pareto margins
#' x <- apply(data, 2, function(i) n/(n + 0.5 - rank(i)))
#' ## Define indices in which we evaluate the estimator
#' indices <- selectGrid(cst = c(0,0.5,1), d = 3)
#' EstimationMaxLinear(x, indices, k = 100, method = "WLS", startingValue = c(0.3,0.5,0.9))
#' indices <- selectGrid(cst = c(0,1), d = 3)
#' EstimationMaxLinear(x, indices, k = 100, method = "Mestimator", startingValue = c(0.3,0.5,0.9))
EstimationMaxLinear <- function(x, indices, k, method, Bmatrix = NULL, Ldot = NULL, biascorr = FALSE, k1 = (nrow(x) - 10), tau = 5,
                         startingValue, Omega = diag(nrow(indices)), iterate = FALSE, covMat=TRUE, GoFtest = FALSE, dist = 0.01, EURO = FALSE) {
    if(EURO){
        if(!is.null(Bmatrix)){
            warning('If EURO = TRUE, Bmatrix is set automatically; the Bmatrix passed by the user is ignored')
        }
        Bmatrix <- BmatrixEURO
    }
    if(is.null(Bmatrix)){
        Bmatrix <- function(pars){return(cbind(pars,1-pars))}
    }
    Bpars <- Bmatrix(startingValue)
    if(any(rowSums(Bpars) > 1) || any(c(Bpars) < 0)){
        warning("invalid parameter values")
    } else if(method == "Mestimator" && biascorr == TRUE){
        warning("biascorr can only be used for method = WLS")
    } else{
        if(method == "Mestimator"){
            return(MestimatorML(x, indices, k, Bmatrix = Bmatrix ,startingValue = startingValue))
        } else if(method == "WLS"){
            if(EURO){
                if(!is.null(Ldot)){
                    warning('If EURO = TRUE, Ldot is set automatically; the Ldot passed by the user is ignored')
                }
                Ldot <- LdotEURO
            }
            if(is.null(Ldot)){
                Ldot <- function(indices,pars){return(phidotBasic(indices,pars))}
            }
            return(WLSestimatorML(x, indices, k, Bmatrix = Bmatrix, phidot = Ldot,
                                  biascorr = biascorr, k1 = k1, tau = tau,
                        startingValue = startingValue,Omega = Omega,iterate=iterate,covMat=covMat,
                        GoFtest = GoFtest, dist = dist))
        }
    }
}


