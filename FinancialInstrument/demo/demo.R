require(FinancialInstrument)

# currencies need to be defined first

currency("USD")
currency("GBP")
currency("EUR")
currency("JPY")

# now exchange rates
exchange_rate("GBPUSD","USD","GBP")
exchange_rate("EURUSD","USD","EUR")
exchange_rate("EURGBP","GBP","EUR")
exchange_rate("EURJPY","JPY","EUR")
exchange_rate("USDJPY","JPY","USD")

# now some stocks 
stock("IBM","USD")
stock("SPY","USD")
stock("DIA","USD")

#Contract specs for options on IBM
option(".IBM","USD",multiplier=100,tick_size=.01, underlying_id="IBM")

#Define tradeable option_series instrument 
# (the expiration date in the suffix_id is a Saturday. 
# You may want to use the Friday before that in 'expires')
option_series(root_id="IBM", suffix_id="110716C175", expires="2011-07-15", callput='call',underlying_id='IBM')
# Also note that even though we defined the root with a primary_id of '.IBM', 
# 'option_series' will still be able to find it when we call it with root_id='IBM' 
option_series(primary_id='IBM_110917P175', expires='2011-09-16')
option_series('IBM_110917C175') #magically figures everthing out...however, expires will be '2011-09' with no day.

#Or use yahoo to help define the specs and all near-term options
option_series.yahoo('SPY')
#option_series.yahoo("SPY",Exp=NULL) # would define all options on SPY

#load.instruments(system.file("data/currencies.csv",package='FinancialInstrument'))
#load.instruments(system.file("data/root_contracts.csv",package='FinancialInstrument'))
#load.instruments(system.file("data/future_series.csv",package='FinancialInstrument'))

# Define a futures root
future("ES",'USD',multiplier=50)
future_series(root_id='ES', suffix_id='Z11', expires='2011-12-16')
# or use magic if you are okay with an 'expires' that is just YYYY-MM
future_series("ES_U11")
# future_series("ESM2") #this would work, but it is recommended to use the underscore to avoid ambiguity.
future_series("ES_M2")
# You can later update the expiration date manually if you like
instrument_attr("ES_M2", 'expires', '2012-06-15')

# bond & bond future

# non-US
stock("BMW","EUR")
BMW <- getSymbols("BMW.DE",src='yahoo',auto.assign=FALSE)
EURUSD <- getSymbols("DEXUSEU",src='FRED',auto.assign=FALSE)

BMW.USD <- redenominate("BMW") #convert prices from EUR to USD 
# Define a synthetic instrument that is BMW denominated in USD instead of EUR.
synthetic("BMW.USD","USD",1,members=c("BMW","EURUSD"))

# Define a spread
getSymbols(c("SPY","DIA")) #download the data for both legs
SPYDIA.fSB <- fn_SpreadBuilder("SPY","DIA", auto.assign=FALSE) #build a 2 leg spread with multiple columns
#or define the spread first
spread("SPY.DIA", "USD", members=c("SPY","DIA"), memberratio=c(1,-1))
SPYDIA.bS <- buildSpread("SPY.DIA", auto.assign=FALSE) #and build it (could be multiple-leg)

SPYDIA.rat <- buildRatio(c("SPY","DIA")) #calculate ratios of prices


##Look at what has been defined.
#bottom-up
buildHierarchy( c("IBM","SPY",".IBM"), c("currency", "multiplier"))
#top-down
it <- instrument.table()
head(it)

it2 <- instrument.table( ,attrs.of='USD') #only show attributes that instrument "USD" also has
head(it2)

