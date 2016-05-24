#' Calculates the Exposure at Default for the IRDs + Commodity example as given in the Basel III regulatory paper
#' @title IRDs+Commodity Example
#' @return The exposure at default (expected value  based on the Basel paper is 936)
#' @export
#' @author Tasos Grivas <tasos@@openriskcalculator.com>
#' @references Basel Committee: The standardised approach for measuring counterparty credit risk exposures
#' http://www.bis.org/publ/bcbs279.htm

ExampleIRDCredit =function()
{

tr1 = CreditSingle(Notional=10000,MtM=20,Currency="USD",Si=0,Ei=3,BuySell='Buy',SubClass='AA',RefEntity='FirmA')
tr2 = CreditSingle(Notional=10000,MtM=-40,Currency="EUR",Si=0,Ei=6,BuySell='Sell',SubClass='BBB',RefEntity='FirmB')
tr3 = CreditIndex(Notional=10000,MtM=0,Currency="USD",Si=0,Ei=5,BuySell='Buy',SubClass='IG',RefEntity='CDX.IG')
tr4 = IRDSwap(Notional=10000,MtM=30,Currency="USD",Si=0,Ei=10,BuySell='Buy')
tr5 = IRDSwap(Notional=10000,MtM=-20,Currency="USD",Si=0,Ei=4,BuySell='Sell')
tr6 = IRDSwaption(Notional=5000,MtM=50,Currency="EUR",Si=1,Ei=11,BuySell='Sell',OptionType='Put',UnderlyingPrice=0.06,StrikePrice=0.05)

trades= list(tr1,tr2,tr3,tr4,tr5,tr6)


# calculating the Exposure-at-Default
EAD = runExampleCalcs(trades)

return(EAD)
}