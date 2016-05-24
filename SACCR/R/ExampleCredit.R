#' Calculates the Exposure at Default for the Credit example as given in the Basel III regulatory paper
#' @title Credit Products Example
#' @return The exposure at default (expected value based on the Basel paper is 381)
#' @export
#' @author Tasos Grivas <tasos@@openriskcalculator.com>
#' @references Basel Committee: The standardised approach for measuring counterparty credit risk exposures
#' http://www.bis.org/publ/bcbs279.htm

ExampleCredit =function()
{

tr1 = CreditSingle(Notional=10000,MtM=20,Currency="USD",Si=0,Ei=3,BuySell='Buy',SubClass='AA',RefEntity='FirmA')
tr2 = CreditSingle(Notional=10000,MtM=-40,Currency="EUR",Si=0,Ei=6,BuySell='Sell',SubClass='BBB',RefEntity='FirmB')
tr3 = CreditIndex(Notional=10000,MtM=0,Currency="USD",Si=0,Ei=5,BuySell='Buy',SubClass='IG',RefEntity='CDX.IG')

trades= list(tr1,tr2,tr3)


# calculating the Exposure-at-Default
EAD = runExampleCalcs(trades)

return(EAD)

}