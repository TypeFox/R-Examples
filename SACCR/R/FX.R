######################################################################
# Create the base FX class
#
# This is used to represent the FX asset class and it is the parent of the FXSwap subclasses
#' @include Trade.R
FX = setRefClass("FX",
                   fields = list(ccyPair = "character"),
                   contains="Trade",
                   methods = list(
                     initialize = function(...){
                       SubClass<<-' '
                       callSuper(...,TradeGroup='FX')
                     }
                   ))
#' 
#' Creates an FX Swap object with the relevant info needed to calculate the Exposure-at-Default (EAD)
#' @title FX Swap Class
#' @param Notional The notional amount of the trade
#' @docType NULL
#' @param MTM      The mark-to-market valuation of the trade
#' @param Currency The currency set that the trade belongs to
#' @param Si The number of years that the trade will take to start (zero if already started)
#' @param Ei The number of years that the trade will expire
#' @param BuySell Takes the values of either 'Buy' or 'Sell'
#' @return An object of type FXSwap
#' @export
#' @author Tasos Grivas <tasos@@openriskcalculator.com>
#' @references Basel Committee: The standardised approach for measuring counterparty credit risk exposures
#' http://www.bis.org/publ/bcbs279.htm
#' @examples
#' 
#'
#' tr1 = FXSwap(Notional=10000,MtM=30,ccyPair="EUR/USD",Si=0,Ei=10,BuySell='Buy')

FXSwap = setRefClass("FXSwap",
                       
                       contains="FX",
                       
                       methods = list(
                         initialize = function(...){
                           callSuper(...,TradeType='Swap')
                         }
                        
                       ))

