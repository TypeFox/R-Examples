#' @title Variable selection in Partial Least Squares
#' 
#' @description A large collection of variable selection methods for use with
#' Partial Least Squares. These include all methods in Mehmood et al. 2012
#' and more.
#' 
#' @importFrom stats approx chisq.test drop.terms fitted formula
#'  median model.frame model.matrix model.response napredict
#'  optim optimise optimize pf predict pt qf qt quantile
#'  residuals runif sd t.test update var wilcox.test
#'  
#' @seealso \code{\link{VIP}} (SR/sMC/LW/RC), \code{\link{filterPLSR}}, \code{\link{spa_pls}}, 
#' \code{\link{stpls}}, \code{\link{truncation}}, \code{\link{bve_pls}}, \code{\link{mcuve_pls}},
#' \code{\link{ipw_pls}}, \code{\link{ga_pls}}, \code{\link{rep_pls}}.
#' 
#' @references T. Mehmood, K.H. Liland, L. Snipen, S. Sæbø, A review of variable selection 
#' methods in Partial Least Squares Regression, Chemometrics and Intelligent Laboratory Systems
#' 118 (2012) 62-69.
#' 
#' @name plsVarSel
NULL
