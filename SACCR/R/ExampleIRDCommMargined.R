#' Calculates the Exposure at Default for the margined IRDs + Commodity example as given in the Basel III regulatory paper
#' @title Margined IRDs+Commodity Example
#' @return The exposure at default (expected value  based on the Basel paper is 1879)
#' @export
#' @author Tasos Grivas <tasos@@openriskcalculator.com>
#' @references Basel Committee: The standardised approach for measuring counterparty credit risk exposures
#' http://www.bis.org/publ/bcbs279.htm

ExampleIRDCommMargined =function()
{

tr1 = Commodity(Notional=10000,MtM= -50,Si=0,Ei=0.75,BuySell='Buy',SubClass='Energy',commodity_type='Oil/Gas')
tr2 = Commodity(Notional=20000,MtM= -30,Si=0,Ei=2,BuySell='Sell',SubClass='Energy',commodity_type='Oil/Gas')
tr3 = Commodity(Notional=10000,MtM= 100,Si=0,Ei=5,BuySell='Buy',SubClass='Metals',commodity_type='Silver')
tr4 = IRDSwap(Notional=10000,MtM=30,Currency="USD",Si=0,Ei=10,BuySell='Buy')
tr5 = IRDSwap(Notional=10000,MtM=-20,Currency="USD",Si=0,Ei=4,BuySell='Sell')
tr6 = IRDSwaption(Notional=5000,MtM=50,Currency="EUR",Si=1,Ei=11,BuySell='Sell',OptionType='Put',UnderlyingPrice=0.06,StrikePrice=0.05)

trades= list(tr1,tr2,tr3,tr4,tr5,tr6)

coll = CSA(thres_cpty = 0, MTA_cpty = 5, IM_cpty = 150, remargin_freq = 5)

MF = coll$CalcMF()

current_collateral = 200

# calculating the add-on
Addon_Aggregate <- CalcAddon(trades, MF)
# calculating the RC and the V-c amount
rc_values <- CalcRC(trades, coll, current_collateral)
# calculating the PFE after multiplying the addon with a factor if V-C<0
PFE <- CalcPFE(rc_values$V_C, Addon_Aggregate)
# calculating the Exposure-at-Default
EAD <- CalcEAD(rc_values$RC,PFE)


return(EAD)
}