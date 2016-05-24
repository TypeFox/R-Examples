#' Calculates the Exposure at Default for the Commodities example as given in the Basel III regulatory paper
#' @title Commodities Example
#' @return The exposure at default (expected value based on the Basel paper is 5406)
#' @export
#' @author Tasos Grivas <tasos@@openriskcalculator.com>
#' @references Basel Committee: The standardised approach for measuring counterparty credit risk exposures
#' http://www.bis.org/publ/bcbs279.htm
 
ExampleComm =function()
{

tr1 = Commodity(Notional=10000,MtM= -50,Si=0,Ei=0.75,BuySell='Buy',SubClass='Energy',commodity_type='Oil/Gas')
tr2 = Commodity(Notional=20000,MtM= -30,Si=0,Ei=2,BuySell='Sell',SubClass='Energy',commodity_type='Oil/Gas')
tr3 = Commodity(Notional=10000,MtM= 100,Si=0,Ei=5,BuySell='Buy',SubClass='Metals',commodity_type='Silver')

trades= list(tr1,tr2,tr3)

# calculating the Exposure-at-Default
EAD = runExampleCalcs(trades)

return(EAD)
}