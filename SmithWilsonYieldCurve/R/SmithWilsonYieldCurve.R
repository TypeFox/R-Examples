#' A package to fit yield curves using the Smith-Wilson method
#' 
#' The main function exposed in this package is fFitSmithWilsonYieldCurve, which takes market data in the form of a vector of cashflow times, a matrix of cashflows and a vector of market prices. It returns an object of class "SmithWilsonYieldCurve". 
#' 
#' A convenience function fFitSmithWilsonYieldCurveToInstruments takes a dataframe containing market instrument data as type, tenor, frequency and rate. It extracts the required vectors and matrices and then calls fFitSmithWilsonYieldCurve.
#' 
#' Objects of class SmithWilsonYieldCurve are a list, the first element of which is a function P(t), which returns the zero coupon bond price of the fitted curve at time t.
#' 
#' Details including mathematics at \url{http://www.not-normal-consulting.co.uk}, or check the EIOPA document in references.
#' 
#' @name SmithWilsonYieldCurve-package
#' @aliases SmithWilsonYieldCurve
#' @docType package
#' @title Fit yield curves using the Smith-Wilson method
#' @author Phil Joubert \email{phil.joubert@@not-normal-consulting.co.uk}
#' @references \url{http://eiopa.europa.eu/fileadmin/tx_dam/files/consultations/QIS/QIS5/ceiops-paper-extrapolation-risk-free-rates_en-20100802.pdf}
#' @examples
#' 	dfInstruments <- data.frame(c("SWAP", "SWAP"), c(1,10), c(1,1), c(0.025, 0.05))
#' 	colnames( dfInstruments ) <- c( "Type", "Tenor", "Frequency", "Rate" )
#' 	Curve <- fFitSmithWilsonYieldCurveToInstruments( dfInstruments, 0.04, 0.1 )
#' 	plot( Curve )
NA