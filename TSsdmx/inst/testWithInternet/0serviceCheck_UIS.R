######################## Unesco Institute for Statistics #######################

require("RJSDMX")

#  sdmxHelp()

#  getFlows('UIS') # 6  eg DEMO_DS   Demogr & socio-econ
#  NY_GDP_PCAP_CD: GDP per capita (current US$)
  z <- getSDMX('UIS', 'DEMO_DS.NY_GDP_PCAP_CD.CAN') # BUG message freq null but should be annual
  if(start(z[[1]]) !=  1970)  stop("test 1 start date changed.")
  if(frequency(z[[1]]) !=  1) stop("test 1  frequency changed.")

  z <- getSDMX('UIS', 'DEMO_DS.NY_GDP_PCAP_CD.CAN+USA+MEX') 

  if(start(z[[1]]) !=  1970)  stop("test 2 start date changed.")
  if(frequency(z[[1]]) !=  1) stop("test 2  frequency changed.")
  if(length(z) !=  3)  stop("test 2  number of series changed.")
