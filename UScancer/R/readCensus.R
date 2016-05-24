readCensus <- function(dsn="/",layer="County_2010Census_DP1",state,fips,codes) {
  usPopFull = readOGR(dsn,layer)
  #data(dplookup)
    
  # create usPop with only the columns and rows we need
  usPop = usPopFull
  
  usPop@data = usPop@data[,as.character(codes$code)]
  names(usPop@data) = codes$desc
  usPop$fips = as.character(usPop$fips)
  usPop@data <- cbind(usPop@data,fips[usPop$fips,c("state","county")])

  # select only those rows that match 'state'
  if (missing(state)) {
    print("Note: state not supplied, including all states")
  } else {
    usPop = usPop[(tolower(usPop$state) %in% tolower(state)),]
  }
  
  return(usPop)
}


# Test
#load("../dplookup.Rdata")
#dsn="../../../CCO/usCensus"
#layer="County_2010Census_DP1"
#usPop <- readCensus(dsn,layer,dplookup)




# Read census subdivision
# data(dplookup)
# data(fipslookup)
# dsn="../../../CCO/usCensus/CouSub_2010Census_DP1"
# layer="CouSub_2010Census_DP1"
# state="Kentucky"
# fips=fipslookup
# codes=dplookup
