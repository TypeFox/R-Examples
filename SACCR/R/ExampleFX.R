#' Calculates the Exposure at Default for the FX product type
#' @title FX Example
#' @return The exposure at default
#' @export
#' @author Tasos Grivas <tasos@@openriskcalculator.com>
#' @references Basel Committee: The standardised approach for measuring counterparty credit risk exposures
#' http://www.bis.org/publ/bcbs279.htm

ExampleFX =function()
{

tr1 = FXSwap(Notional=10000,MtM=30,ccyPair="EUR/USD",Si=0,Ei=10,BuySell='Buy')
tr2 = FXSwap(Notional=20000,MtM=-20,ccyPair="EUR/USD",Si=0,Ei=4,BuySell='Sell')
tr3 = FXSwap(Notional=5000,MtM=50,ccyPair="GBP/USD",Si=1,Ei=11,BuySell='Sell')

trades= list(tr1,tr2,tr3)

# calculating the Exposure-at-Default
EAD = runExampleCalcs(trades)

return(EAD)

}