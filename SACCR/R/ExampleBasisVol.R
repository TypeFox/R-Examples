#' Calculates the Exposure at Default for a trade set containing basis and volatility transactions
#' @title Basis+Volatility trades Example
#' @return The exposure at default
#' @export
#' @author Tasos Grivas <tasos@@openriskcalculator.com>
#' @references Basel Committee: The standardised approach for measuring counterparty credit risk exposures
#' http://www.bis.org/publ/bcbs279.htm

ExampleBasisVol =function()
{

# creating the trade objects and storing them in a list
tr1 = IRDSwap(Notional=10000,MtM=30,Currency="USD",Si=0,Ei=10,BuySell='Buy',pay_leg_type='Float', pay_leg_ref  = "CDOR",
              pay_leg_tenor= "1M",rec_leg_type='Float',rec_leg_ref="CORRA",rec_leg_tenor="3M")
tr2 = CommSwap(Notional=10000,MtM=-20,Currency="USD",Si=0,Ei=4,BuySell='Sell',SubClass='Energy',commodity_type='Oil/Gas',pay_leg_type='Commodity', pay_leg_ref  = "Brent", rec_leg_type='Commodity', rec_leg_ref  = "Gas")
tr3 = IRDSwapVol(Notional=5000,MtM=50,Currency="EUR",Si=1,Ei=11,BuySell='Sell', reference  = "CDOR", vol_strike=0.2, annualization_factor=252)
tr4 = IRDSwap(Notional=10000,MtM=30,Currency="USD",Si=0,Ei=10,BuySell='Buy')

trades= list(tr1,tr2,tr3,tr4)

# calculating the Exposure-at-Default
EAD = runExampleCalcs(trades)

return(EAD)
}