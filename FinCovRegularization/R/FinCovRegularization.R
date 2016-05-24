#' FinCovRegularization: Covariance Matrix Estimation and Regularization 
#'   for Finance
#'
#' Estimation and regularization for covariance matrix of asset returns. 
#' For covariance matrix estimation, three major types of factor models are 
#' included: macroeconomic factor model, fundamental factor model and 
#' statistical factor model. For covariance matrix regularization, 
#' four regularized estimators are included: banding, tapering, 
#' hard-thresholding and soft-thresholding. The tuning parameters of these 
#' regularized estimators are selected via cross-validation.
#'
#' @docType package
#' @name FinCovRegularization
NULL

#' 10 stock and S&P 500 excess returns
#'
#' A dataset containing monthly excess returns of 10 stocks and 
#' S$P 500 index return from January 1990 to December 2003
#'
#' @docType data
#' @format A matrix with 168 rows and 11 variables
#' @usage data(m.excess.c10sp9003)
"m.excess.c10sp9003"