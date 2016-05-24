########### URLS for Canada scraper
 
STATION.URL  <- "https://scraperwiki.com/api/1.0/datastore/sqlite?format=csv&name=can-weather-stations&query=select+*+from+`swdata`"
BASE.URL <-"http://climate.weatheroffice.gc.ca/climateData/bulkdata_e.html?timeframe=3&Prov=XX&StationID="
YEAR.URL <-"&Year="
FORMAT.URL <-"&Month=1&Day=1&format=csv"
MASTER.STATION.LIST <- "EnvCanadaMaster.csv"
MONTHLY.STATION.LIST <- "MonthlyStations.Env.csv"
STARTING.STATION.ID  <- 99111111
DATA.DIRECTORY       <- "DataDirectory"