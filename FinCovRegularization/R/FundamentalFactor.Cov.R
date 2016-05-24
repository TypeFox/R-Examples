#' @title Covariance Matrix Estimation by Fundamental Factor Model
#'
#' @description
#' Estimate covariance matrix by fitting a fundamental factor model
#' using OLS or WLS regression
#'
#' @param assets a N*p matrix of asset returns, N indicates sample size
#'   and p indicates the dimension of asset returns
#' @param exposure a p*q matrix of exposure indicator for the
#'   fundamental factor model, p corresponds to the dimension of
#'   asset returns, q indicates the number of fundamental industries
#' @param method a character, indicating regression method: "OLS" or "WLS"
#' @return an estimated p*p covariance matrix
#' @examples
#' data(m.excess.c10sp9003)
#' assets <- m.excess.c10sp9003[,1:10]
#' Indicator <- matrix(0,10,3)
#' dimnames(Indicator) <- list(colnames(assets),c("Drug","Auto","Oil"))
#' Indicator[c("ABT","LLY","MRK","PFE"),"Drug"] <- 1
#' Indicator[c("F","GM"),"Auto"] <- 1
#' Indicator[c("BP","CVX","RD","XOM"),"Oil"] <- 1
#' FundamentalFactor.Cov(assets,exposure=Indicator,method="WLS")
#' @importFrom stats cov
#' @importFrom stats var
#' @export

FundamentalFactor.Cov <- function(assets, exposure, method = "WLS") {
  if ((method %in% c("OLS","WLS")) == FALSE) {
    stop("This function only support two methods: OLS and WLS")
  }
  # Multivariate OLS regression to estimate OLS factor returns
  F.hat.ols <- solve(crossprod(exposure)) %*% t(exposure) %*% t(assets)
  # Compute matrix of OLS industry factor model residuals
  E.hat.ols <- assets - t(exposure %*% F.hat.ols)
  # Compute OLS residual variances from time series of errors
  diagD.hat.ols  <-  apply(E.hat.ols, 2, var)
  if (method == "OLS") {
    # Compute sample covariance matrix of OLS estimated factors
    COV.OLS <- exposure %*% cov(t(F.hat.ols)) %*% t(exposure) + diag(diagD.hat.ols)
    dimnames(COV.OLS) <- list(colnames(assets), colnames(assets))
    return(COV.OLS)
  } else if (method == "WLS") {
    # Multivariate WLS regression to estimate WLS factor returns
    Dinv.hat.ols <- diag(diagD.hat.ols^(-1))
    # Compute factor mimicking portfolios weights matrix
    H.hat <- solve(t(exposure) %*% Dinv.hat.ols %*% exposure) %*% t(exposure) %*% Dinv.hat.ols
    F.hat.wls <- as.matrix(assets) %*% t(H.hat)
    # Compute matrix of WLS industry factor model residuals
    E.hat.wls <- assets - t(exposure %*% t(F.hat.wls))
    # Compute OLS residual variances from time series of errors
    diagD.hat.wls  <-  apply(E.hat.wls, 2, var)
    # Compute sample covariance matrix of WLS estimated factors
    COV.WLS <- exposure %*% cov(F.hat.wls) %*% t(exposure) + diag(diagD.hat.wls)
    dimnames(COV.WLS) <- list(colnames(assets), colnames(assets))
    return(COV.WLS)
  }
}
