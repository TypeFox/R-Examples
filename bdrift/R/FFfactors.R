#' @title Dataset to Estimate Multi-Factor Models for Return Samples
#'
#' @description A dataset containing simple daily returns of ExxonMobile,
#' BlackRock's Large Cap Core Inv A fund, and Vanguard 500 ETF fund
#' as well as all necessary factor data from Kenneth French's data library to estimate
#' a five-factor model. The dataset contains data from Sep-09-2010 until Nov-30-2015.
#' 
#' @docType data
#' @keywords datasets
#' @name FFfactors
#' @usage FFfactors
#' @format A xts object with 1316 rows and 10 variables. The variables
#' are as follows:
#' \describe{
#' \item{XOM}{simple daily returns of ExxonMobile (NYSE:XOM) less the risk-free rate}
#' \item{MDLRX}{simple daily returns of BlackRock Large Cap Core Inv A fund (MDLRX)
#'  less the risk-free rate}
#' \item{VOO}{simple daily returns of Vanguard 500 ETF (VOO) less the risk-free rate}
#' \item{SP500}{simple daily returns of Standard&Poor's 500 index less the risk-free rate}
#' \item{Mkt.RF}{simple daily returns of the US market (all NYSE, AMEX, and NASDAQ
#' firms) less the risk-free rate}
#' \item{SMB}{daily small-minus-big factor(size factor)}
#' \item{HML}{daily high-minus-low factor(value factor)}
#' \item{RMW}{daily robust-minus-weak factor(profitability factor)}
#' \item{CMA}{daily conservative-minus-aggressive factor(investment factor)}
#' \item{RF}{daily risk-free rate}
#' }
#' @source \code{XOM}, \code{MDLRX}, \code{VOO}, and \code{SP500} was
#' retrieved rom Yahoo Finance via \code{\link[quantmod]{getSymbols}}
#' from the \code{\link[quantmod]{quantmod}} package.
#' All other factors were retrieved from Kenneth French's data library
#' \emph{Kenneth French's data library} at \url{
#' http://mba.tuck.dartmouth.edu/pages/faculty/ken.french/data_library.html}.
#' \code{RF} was originally provided by Ibbotson Associates.
#'
#' @references Fama, E. F., & French, K. R. (1993).
#' Common risk factors in the returns on stocks and bonds.
#' \emph{Journal of financial economics}, 33(1), 3-56.
#'
#' Fama, E. F., & French, K. R. (2015).
#' A five-factor asset pricing model.
#' \emph{Journal of financial economics}, 116(1), 1-22.
NULL
