#' Creates a Commodity Object with the relevant info needed to calculate the Exposure-at-Default (EAD)
#' @title Commodity Class
#' @include Trade.R
#' @param Notional The notional amount of the trade
#' @param MTM      The mark-to-market valuation of the trade
#' @param Currency The currency set that the trade belongs to
#' @param Si The number of years that the trade will take to start (zero if already started)
#' @param Ei The number of years that the trade will expire
#' @param BuySell Takes the values of either 'Buy' or 'Sell'
#' @param commodity_type Takes the values of 'Oil/Gas','Silver','Electricity' etc.
#' @return An object of type Commodity
#' @export
#' @author Tasos Grivas <tasos@@openriskcalculator.com>
#' @references Basel Committee: The standardised approach for measuring counterparty credit risk exposures
#' http://www.bis.org/publ/bcbs279.htm
#' @examples
#'
#' ## the Commodity trade given in the Basel regulation Credit example
#' tr1 = Commodity(Notional=10000,MtM=-50,Si=0,Ei=0.75,
#' BuySell='Buy',SubClass='Energy',commodity_type='Oil/Gas')

Commodity = setRefClass("Commodity",
                         fields = list(commodity_type      = 'character'
                         ),
                      contains="Trade",
                      methods = list(
                        initialize = function(...){
                          callSuper(...,TradeGroup='Commodity')}
                      ))
#' Creates a Commodity Swap Object with the relevant info needed to calculate the Exposure-at-Default (EAD)
#' @title Commodity Swap Class
#' @include Swap.R
#' @return An object of type CommSwap
#' @export
#' @author Tasos Grivas <tasos@@openriskcalculator.com>
#' @references Basel Committee: The standardised approach for measuring counterparty credit risk exposures
#' http://www.bis.org/publ/bcbs279.htm
CommSwap = setRefClass("CommSwap",
                       fields = list(),
                       contains=c("Swap","Commodity"),
                       methods = list(
                         initialize = function(...){
                           callSuper(...)}
                       ))
