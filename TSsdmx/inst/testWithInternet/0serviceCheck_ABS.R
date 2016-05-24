########################## Australian Bureau of Statistics ####################
# http://www.abs.gov.au

require("RJSDMX")

# getFlows('ABS')

#  sdmxHelp()

##################  quarterly #################

tts <-  getSDMX('ABS', 'BOP.1.100.10.Q') 
names(tts)

if ("1959 Q3" != start(tts[[1]])) stop("ABS test 1 start date changed.")
if   (4 != frequency(tts[[1]]))   stop("ABS test 1 frequency failure.")

tts <-  getSDMX('ABS', 'CPI.1.50.10001.10.Q', start="1960-Q1", end="2010-Q4") 
names(tts)

if ("1960 Q1" != start(tts[[1]])) stop("ABS test 2 start date failure.")
if ("2010 Q4" != end(tts[[1]]))   stop("ABS test 2  end  date failure.")

#tts <-  getSDMX('ABS', 'CPI.1.*.10001.10.Q') #19 provider BUG (White spaces...)
## tts <-  getSDMX('ABS', 'CPI.1*') 1461 series


#if ("1959 Q3" != start(tts[[1]])) stop("ABS test 1 start date changed.")
#if   (4 != frequency(tts[[1]]))   stop("ABS test 1 frequency failure.")

#"ABS,PPI_1 ; Producer Price Indexes by Industry"
#"ABS,LF ; Labour Force"
#"ABS,BOP ; Balance of Payments"

##################  monthly #################

#tts <-  getSDMX('ABS', 'RT.0*M') # empty
#tts <-  getSDMX('ABS', 'RT.0*')   #178 series

tts <-  getSDMX('ABS', 'RT.0.2.15.10.M', start="2008-05", end="2014-07") 
names(tts)

if ("May 2008" != start(tts[[1]])) stop("ABS test 2 start date failure.")
if ("Jul 2014" != end(tts[[1]]))   stop("ABS test 2  end  date failure.")
if   (12 != frequency(tts[[1]]))   stop("ABS test 2 frequency failure.")
