require(FinancialInstrument)
require(quantmod)

#Similar to the other demo, but demonstrates new functionality

# currencies need to be defined first
#load.instruments(system.file('data/currencies.csv',package='FinancialInstrument')))
currency(c("USD","GBP","EUR","JPY"))

# now exchange rates
exchange_rate(counter_currency=c("GBP","EUR"), currency=c("GBP","JPY","USD"))
exchange_rate("USDJPY")

# now some stocks 
stock(c("IBM","SPY","DIA","AAPL"),"USD")

#Contract specs for options
option(".IBM","USD",multiplier=100,tick_size=0.01, underlying_id="IBM")
#If we don't provide the currency, it will use the currency of the underlying stock
option(".SPY",multiplier=100,tick_size=0.01,underlying_id='SPY')
option(".GOOG",multiplier=100,underlying_id=stock("GOOG","USD")) #underlying can be defined on-the-fly
#if you don't provide a primary_id, it will become the underlying_id with a dot in front of it
option(multiplier=100,underlying_id=stock("GS","USD"))
#don't remember what primary_id was used? Pass the underlying_id to getInstrument and tell it to find an option
getInstrument('GS',type='option')

#Define tradeable option_series instrument 
# (the expiration date in the suffix_id is a Saturday. 
# You may want to use the Friday before that in 'expires')
option_series(root_id="IBM", suffix_id="110716C175", expires="2011-07-15", callput='call',underlying_id='IBM')
# Also note that even though we defined the root with a primary_id of '.IBM', 
# 'option_series' will still be able to find it when we call it with root_id='IBM' 
option_series('IBM_110917C175') #magically figures everthing out...however, since we didn't provide expires, it will be '2011-09' with no day.
option_series(primary_id='IBM_110917P175', expires='2011-09-16')

#Or use yahoo to help define the specs and all near-term options
option_series.yahoo('SPY')
#option_series.yahoo("SPY",Exp=NULL) # would define all options on SPY

#load.instruments(system.file("data/currencies.csv",package='FinancialInstrument'))
#load.instruments(system.file("data/root_contracts.csv",package='FinancialInstrument'))
#load.instruments(system.file("data/future_series.csv",package='FinancialInstrument'))

# Define a futures root (don't have to provide currency if the underlying is defined.)
future("ES",multiplier=50,underlying_id=synthetic("SPX","USD",exchange='CBOE')) #underlying can be defined on the fly.
future_series(root_id='ES', suffix_id='Z11', expires='2011-12-16')
# or use magic if you are okay with an 'expires' that is just YYYY-MM
future_series("ES_U11") 
#or provide expires and let magic take care of the rest
future_series("ES_M2", expires='2012-06-15') 

# You can later update the expiration date manually if you like
getInstrument("ES_U11")
instrument_attr("ES_U11", 'expires', '2011-09-16')
getInstrument("ES_U11")

future_series("ESM2") #this works, but it is recommended to use the underscore to avoid ambiguity.


##Using the src arg.
#if you specify 'src' when defining an instrument, a call will be made to setSymbolLookup 
#so that getSymbols will know where to get the data.
stock("YHOO",'USD',src='google')
synthetic("AORD",currency("AUD"),src=list(name='^AORD',src='yahoo')) #we also are defining "AUD" on-the-fly
exchange_rate("AUDUSD", src=list(name='AUD/USD', src='oanda'))
#stock("GM","USD",src=list(src='FI',dir='/mnt/data') 
#now getSymbols will look in appropriate places for data
getSymbols(c("YHOO","AORD","AUDUSD"), from=Sys.Date()-250)
showSymbols()

# Define a spread
getSymbols(c("SPY","DIA")) #download the data for both legs
SPYDIA.fSB <- fn_SpreadBuilder("SPY","DIA", auto.assign=FALSE) #build a 2 leg spread with multiple columns
fn_SpreadBuilder(stock(c("CVX","XOM"),"USD")) #define stocks, download data, calculate spread, and define and assign spread

#or define the spread first
spread("SPY.DIA", "USD", members=c("SPY","DIA"), memberratio=c(1,-1))

#let 'spread' make the primary_id for you. Since both members have the same currency, you don't have to provide it. 
spread(members=c("AAPL","GOOG"),memberratio=c(1,-1)) 
SPYDIA.bS <- buildSpread("SPY.DIA", auto.assign=FALSE) # (could be multiple-leg.  See VX fly below)
buildSpread('AAPL.GOOG')

#VX fly -- need qmao package for the getSymbols.cfe method
future("VX", "USD", 1000, src='cfe', underlying_id=synthetic("VIX","USD")) #first the root
butterfly(members=future_series(c('VX_F11','VX_G11','VX_H11'),src='cfe')) #we can define the future_series on-the-fly
if (require(qmao)) {
    buildSpread("VX_F11.G11.H11")
    head(VX_F11.G11.H11)
}

#ratios of prices
SPYDIA.rat <- buildRatio(c("SPY","DIA")) #calculate ratios of prices


