######################## Unesco Institute for Statistics #######################

require("TSsdmx")

  uis <- TSconnect("sdmx", dbname="UIS")

  z <- TSget('DEMO_DS.NY_GDP_PCAP_CD.CAN', uis) 
  if(start(z) !=  1970)  stop("test 1 start date changed.")
  if(frequency(z) !=  1) stop("test 1  frequency changed.")

  z <- TSget('DEMO_DS.NY_GDP_PCAP_CD.CAN+USA+MEX', uis) 

  if(start(z) !=  1970)  stop("test 2 start date changed.")
  if(frequency(z) !=  1) stop("test 2  frequency changed.")
  if(tframe::nseries(z) !=  3)  stop("test 2  number of series changed.")
