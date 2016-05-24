require("RJSDMX")


############################ EUROSTAT ############################

#[ http://epp.eurostat.ec.europa.eu/portal/page/portal/eurostat/home ]

#http://epp.eurostat.ec.europa.eu/portal/page/portal/statistics/search_database
#   >Economy and finance
#      >National accounts (including GDP) (ESA95) (na
#         >Quarterly national accounts (namq) 
#              >GDP and main components (namq_gdp)

  names(getDimensions('EUROSTAT','ei_nama_q')) 
  getCodes('EUROSTAT','ei_nama_q', 'FREQ')
  
  nm <- getFlows('EUROSTAT')
  length(nm)  # 5717 on 7 Nov 2014
  
  getFlows('EUROSTAT', "namq_gdp_c")  # length 1

  getFlows('EUROSTAT', "ei_nama_q")  # length 1


#### quarterly ####
  # as of Sept 2015 next fails if compression is enabled (BUG #76)
  # compression can be disbled in .SdmxClient config file.
  tts1 <- getSDMX('EUROSTAT', "ei_nama_q.Q.MIO-EUR.SWDA.CP.NA-P72.IT") 
  names(tts1)

  if("1980 Q1" != start(tts1[[1]]))
                   stop("start test for EUROSTAT quarterly data failed.")
  if(4 != frequency(tts1[[1]])) 
             stop(  "frequency test for EUROSTAT quarterly data failed.")

  tts2 <- getSDMX('EUROSTAT', "ei_nama_q.Q.MIO-EUR.SWDA.CP.NA-P72.IT",
                  start="1990")[[1]]

  if("1990 Q1" != start(tts2))
        stop("EUROSTAT quarterly start specification 2 failure.")

  tts3 <- getSDMX('EUROSTAT', "ei_nama_q.Q.MIO-EUR.SWDA.CP.NA-P72.IT",
	    start="1990-Q1", end="2012-Q2")[[1]]

  if("1990 Q1" != start(tts3))
        stop("EUROSTAT quarterly start specification 3 failure.")
  if("2012 Q2" != end(tts3))
        stop("EUROSTAT quarterly  end  specification 3 failure.")

  #tts2 = getSDMX('EUROSTAT', 'ei_nama_q.Q.*.*.*.*.IT')   # works
  #tts2 = getSDMX('EUROSTAT', 'ei_nama_q.Q.*.*.CP.*.IT')  # works
  #tts2 = getSDMX('EUROSTAT', 'ei_nama_q.Q.*.NSA.CP.*.IT')  # works
  #tts2 = getSDMX('EUROSTAT', 'ei_nama_q.Q.*.*.CP.*.*.*') NO

  #tts2 = getSDMX('EUROSTAT', 'ei_nama_q.Q.MIO-EUR.NSA.*.*.IT') 
       #  above has 84 series Feb 2015, but may change
  #tts2 = getSDMX('EUROSTAT', 'ei_nama_q.Q.MIO-EUR.NSA.CP.*.IT') #  28 series
  #names(tts2)

  #nm[167]   #                "ei_nama_q.Q.MIO-EUR.NSA.CP.NA-P72.IT"
  #nm[168]   #                "ei_nama_q.Q.MIO-EUR.SWDA.CP.NA-P72.IT"

  # for (i in 1: length(tts2)) print( any(! is.nan(tts2[[i]])))
  # for (i in 1: length(tts2)) print( sum(! is.nan(tts2[[i]])))


  # z <- getSDMX('EUROSTAT', 'ei_nama_q.Q.MIO-EUR.NSA.CLV2000.*.IT')[[1]]

  # if("1980 Q1" != start(z)) stop("EUROSTAT quarterly retrieval start changed.")
  # if(4 != frequency(z)) stop("EUROSTAT quarterly retrieval frequency error.")
  
 
