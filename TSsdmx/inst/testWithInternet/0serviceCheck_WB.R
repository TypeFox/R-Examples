######################## World Bank #######################

# See http://data.worldbank.org/developers

require("RJSDMX")

#  sdmxHelp()

#  getFlows('WB') # $WDI

#  getDimensions('WB', 'WDI') 
#  $FREQ     [1] "WB/CL_FREQ_WDI"
#  $SERIES   [1] "WB/CL_SERIES_WDI"
#  $REF_AREA [1] "WB/CL_REF_AREA_WDI"

# names(getCodes('WB', 'WDI', 'FREQ'))

#BM_GSR_MRCH_CD  Goods imports (BoP, current US$)
#BM_GSR_NFSV_CD  Services imports (BoP, current US$)

#The World Bank (beta) has a slightly unconventional indication for time series. #Hopefully this will change in a future releases.
#Even if the declared structure is FREQ.SERIES.REF_AREA, you have to build
# the queries as REF_AREA.SERIES. for example:

  if(! TSsdmx::verifyQuery('WB', 'WDI.*.*.USA'))
     stop("verifyQuery 1a does not verify. Provider changed something.")

  if(! TSsdmx::verifyQuery('WB', 'WDI.A.*.USA'))
     stop("verifyQuery 1b does not verify. Provider changed something.")

  if(! TSsdmx::verifyQuery('WB', 'WDI.A.SP_POP_TOTL.USA'))
     stop("verifyQuery 1c does not verify. Provider changed something.")


  #z = getSDMX('WB', 'WDI/CHN.SP_POP_TOTL')# fails
  #z = getSDMX('WB', 'WDI.CHN.SP_POP_TOTL')# fails
  z  = getSDMX('WB', 'WDI.CHN.SP_POP_TOTL', start='2000', end='2010')# works
  #z = getSDMX('WB', 'WDI/CHN.SP_POP_TOTL', start='2000', end='2010')# works too
 
  if(start(z[[1]]) !=  2000)  stop("test 1 start date changed.")
  if(end(z[[1]])   !=  2010)  stop("test 1  end  date changed.")
  if(frequency(z[[1]]) !=  1) stop("test 1  frequency changed.")
  if(length(z) !=  1)  stop("test 1  number of series changed.")
 

  #z <- getSDMX('WB', 'WDI.CHN.*', start='2000', end='2010') #fails
  #z <- getSDMX('WB', 'WDI.CHN.', start='2000', end='2010') #fails
  #z <- getSDMX('WB', 'WDI/CHN.*', start='2000', end='2010')  #fails
  #z <- getSDMX('WB', 'WDI/CHN.', start='2000', end='2010')  #fails

  #z <- getSDMX('WB', 'WDI.A.BM_GSR_MRCH_CD.CAN')#fails

  #z <- getSDMX('WB', 'WDI.CAN.A.BM_GSR_MRCH_CD')#fails
  #z <- getSDMX('WB', 'WDI.CAN.A.BM_GSR_MRCH_CD')#fails
  #z <- getSDMX('WB', 'WDI.CAN.BM_GSR_MRCH_CD')#fails

  #z  = getSDMX('WB', 'WDI.*.SP_POP_TOTL', start='2000', end='2010')# fails
  #z  = getSDMX('WB', 'WDI.CAN.SP_POP_TOTL', start='2000', end='2010')# fails

  z  = getSDMX('WB', 'WDI.USA.SP_POP_TOTL', start='1975', end='2011')# works

  if(start(z[[1]]) !=  1975)  stop("test 2 start date changed.")
  if(frequency(z[[1]]) !=  1) stop("test 2  frequency changed.")
  if(end(z[[1]])   !=  2011)  stop("test 2  end  date changed.")
  if(length(z) !=  1)  stop("test 2 number of series  changed.")
