if (FALSE){ 
######################## National Bank Belguim #######################

require("RJSDMX")

#  sdmxHelp()
# National Bank Belgium
#  getFlows('NBB') # many (eg $HICP )
  #z <- getSDMX('NBB', 'HICP.000000..') #works
  #z <- getSDMX('NBB', 'HICP.000000.*.*') #works
  #names(z)
  #z <- getSDMX('NBB', 'HICP.000000.BE.') #works
  z <- getSDMX('NBB', 'HICP.000000.BE.*') #works
  z <- getSDMX('NBB', 'HICP.000000.BE.A') #works
  z <- getSDMX('NBB', 'HICP.000000.BE.M') #works

}
