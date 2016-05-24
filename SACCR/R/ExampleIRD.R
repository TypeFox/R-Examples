#' Calculates the Exposure at Default for the IRD example as given in the Basel III regulatory paper
#' @title IRDs Example
#' @return The exposure at default (expected value based on the Basel paper is 569)
#' @export
#' @author Tasos Grivas <tasos@@openriskcalculator.com>
#' @references Basel Committee: The standardised approach for measuring counterparty credit risk exposures
#' http://www.bis.org/publ/bcbs279.htm

ExampleIRD =function()
{


# creating the trade objects and storing them in a list
tr1 = IRDSwap(Notional=10000,MtM=30,Currency="USD",Si=0,Ei=10,BuySell='Buy')
tr2 = IRDSwap(Notional=10000,MtM=-20,Currency="USD",Si=0,Ei=4,BuySell='Sell')
tr3 = IRDSwaption(Notional=5000,MtM=50,Currency="EUR",Si=1,Ei=11,BuySell='Sell',OptionType='Put',UnderlyingPrice=0.06,StrikePrice=0.05)

trades= list(tr1,tr2,tr3)

# calculating the Exposure-at-Default
EAD = runExampleCalcs(trades)

return(EAD)
}