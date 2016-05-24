##' Basic Finance: NPV/IRR/annuities, bond pricing, Black Scholes 
##' 
##' 
##' This package implements the basic financial analysis functions
##' similar to (but not identical to) what is available in most
##' spreadsheet software.  This includes finding the IRR, NPV and
##' duration of possibly irregularly spaced cash flows and
##' annuities. Bond pricing, YTM and duration calculations are
##' included. Black Scholes option pricing, Greeks and implied
##' volatility are also provided.
##'
##' Important functions include:
##' 
##' \code{\link{npv}}, \code{\link{irr}}, \code{\link{duration}}, 
##' \code{\link{annuity.pv}}, \code{\link{bond.price}}, \code{\link{bond.yield}},
##' \code{\link{GenBS}}, \code{\link{GenBSImplied}}
##' 
##' For more details, see the vignette
##' 
##' @name jrvFinance-package
##' @aliases jrvFinance-package jrvFinance
##' @docType package
##' 
##' @author Prof. Jayanth R. Varma \email{jrvarma@@iimahd.ernet.in}
##' @references
##' The 30/360 day count was converted from C++ code in the QuantLib library.
##' The Newton Raphson solver was converted from C++ code in the Boost library

NULL



