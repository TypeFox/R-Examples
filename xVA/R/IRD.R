######################################################################
# Create the base IRD class
#
# This is used to represent the Interest Rate Derivatives product group and it is the parent of the
# subclasses: IRDSwap
#' @include Trade.R

IRD = setRefClass("IRD",
                   # the timebuckets grouping is only relevant for IRDs
                   fields = list(TimeBucket = "numeric"),
                   contains="Trade",
                   methods = list(
                     initialize = function(...){
                       SubClass<<-' '
                       callSuper(...,TradeGroup='IRD')
                     },
                     SetTimeBucket = function() {
                       ## sets the time bucket based on the maturity of the trade
                       if(Ei<1)         timebucket = 1;
                       if(Ei>=1&&Ei<=5) timebucket = 2;
                       if(Ei>5)         timebucket = 3;
                       return(timebucket)
                     }
                   ))
#' Creates an IR Swap Object with the relevant info needed to calculate the Exposure-at-Default (EAD)
#'
#' @title IR Swap Class
#' @param Notional The notional amount of the trade
#' @param MTM      The mark-to-market valuation of the trade
#' @param Currency The currency set that the trade belongs to
#' @param Si The number of years that the trade will take to start (zero if already started)
#' @param Ei The number of years that the trade will expire
#' @param BuySell Takes the values of either 'Buy' or 'Sell'
#' @param swap_rate The rate of the fixed leg of the swap
#' @return An object of type IRSwap
#' @export
#' @examples
#'
#' # the IR Swap trade given in the Basel regulation IR example
#' tr1 = IRSwap(Notional=10000,MtM=30,Currency="USD",Si=0,Ei=10,BuySell='Buy')

IRSwap = setRefClass("IRSwap",
                      fields = list(swap_rate = 'numeric'),
                       contains="IRD",

                       methods = list(
                         initialize = function(...){
                           callSuper(...,TradeType='Swap')
                         }
                       ))
