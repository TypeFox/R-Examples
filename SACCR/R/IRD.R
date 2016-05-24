######################################################################
# Create the base IRD class
#
# This is used to represent the Interest Rate Derivatives product group and it is the parent of the
# subclasses: IRDSwap and IRDSwaption
#' @include Trade.R

IRD = setRefClass("IRD",
                   # the timebuckets grouping is only relevant for IRDs
                   fields = list(TimeBucket = "numeric"),
                   contains=c("Trade"),
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

#' @title IRD Swaption Class
#' @description  Creates an IRD Swaption Object with the relevant info needed to calculate the Exposure-at-Default (EAD)
#' @param Notional The notional amount of the trade
#' @param MTM      The mark-to-market valuation of the trade
#' @param Currency The currency set that the trade belongs to
#' @param Si The number of years that the trade will take to start (zero if already started)
#' @param Ei The number of years that the trade will expire
#' @param BuySell Takes the values of either 'Buy' or 'Sell'
#' @param OptionType      Takes the values of either 'Put' or 'Call'
#' @param UnderlyingPrice The current price of the underlying
#' @param StrikePrice     The strike price of the option
#' @return An object of type IRDSwaption
#' @export
#' @author Tasos Grivas <tasos@@openriskcalculator.com>
#' @references Basel Committee: The standardised approach for measuring counterparty credit risk exposures
#' http://www.bis.org/publ/bcbs279.htm
#' @examples
#'
#' # the Swaption trade given in the Basel regulation IRD example
#' tr3 = IRDSwaption(Notional=5000,MtM=50,Currency="EUR",Si=1,Ei=11,BuySell='Sell',
#' OptionType='Put',UnderlyingPrice=0.06,StrikePrice=0.05)



IRDSwaption = setRefClass("IRDSwaption",
                       # Swaption contains unique fields relevant to the option's features
                       fields = list(OptionType      = 'character',
                                     UnderlyingPrice = 'numeric',
                                     StrikePrice     = 'numeric'
                       ),
                       contains="IRD",

                       methods = list(
                         initialize = function(...){
                            callSuper(...,TradeType='Option')
                         }
                       ))
#' Creates an IRD Swap Object with the relevant info needed to calculate the Exposure-at-Default (EAD)
#'
#' @title IRD Swap Class
#' @param Notional The notional amount of the trade
#' @param MTM      The mark-to-market valuation of the trade
#' @param Currency The currency set that the trade belongs to
#' @param Si The number of years that the trade will take to start (zero if already started)
#' @param Ei The number of years that the trade will expire
#' @param BuySell Takes the values of either 'Buy' or 'Sell'
#' @return An object of type IRDSwap
#' @export
#' @examples
#'
#' # the IRD Swap trade given in the Basel regulation IRD example
#' tr1 = IRDSwap(Notional=10000,MtM=30,Currency="USD",Si=0,Ei=10,BuySell='Buy')

IRDSwap = setRefClass("IRDSwap",

                       contains=c("IRD","Swap"),

                       methods = list(
                         initialize = function(...){
                           callSuper(...,TradeType='Swap')
                         }
                       ))

#' Creates an IRD Swap Volatility-based Object with the relevant info needed to calculate the Exposure-at-Default (EAD)
#'
#' @title IRD Swap Volatility Class
#' @include Vol.R
#' @include Swap.R
#' @return An object of type IRDSwapVol
#' @export

IRDSwapVol = setRefClass("IRDSwapVol",

                      contains=c("IRD","Swap","Vol"),

                      methods = list(
                        initialize = function(...){
                          callSuper(...,TradeType='Swap')
                        }
                      ))
