########################### "INEGI" ###########################
##### Instituto Nacional de Estadistica y Geografia (Mexico) ######

#For further info about the provider content you may want to check the INEGI SDMX page:
#http://www.inegi.org.mx/inegi/contenidos/servicios/sdmx/

#  sdmxHelp()

require("RJSDMX")

#  getFlows('INEGI') 

getDimensions('INEGI','DF_STES') 
# FREQ   REFERENCE_AREA  SUBJECT  MEASURE  ADJUSTMENT UNIT

getCodes('INEGI','DF_STES', 'FREQ')

nm <- getCodes('INEGI','DF_STES', 'REFERENCE_AREA')
nm$MX
sb <- getCodes('INEGI','DF_STES', 'SUBJECT')
sb$MANMM101 # code exists but series does not
nm <- getCodes('INEGI','DF_STES', 'MEASURE')
nm$ST
nm <- getCodes('INEGI','DF_STES', 'ADJUSTMENT')
nm$NSA
nm <- getCodes('INEGI','DF_STES', 'UNIT')
nm$MXN

# note: all DF_STES series are monthly (as of Feb12,2015)

#tts <- getSDMX("INEGI", 'DF_STES.*.*.*.*.*.*') 
#names(tts)  #26

tts <- getSDMX("INEGI", 'DF_STES.M.MX.*.*.*.*') 
  
if ("Jan 1991" != start(tts[[1]]))  stop("INEGI test 1 start date error.")
if (12 != frequency(tts[[1]]))      stop("INEGI test 1 frequency error.")
if (26 != length(names(tts)))     stop("INEGI test 1 number of series changed.")

names(tts)  #26
sb[c("SLRTCR03", "XTEXVA0", "LMUNRLTT", "ULQCWS02", "LCEAMN04", "HOHWMN03", "PRMNTO01")]


if(! 'DF_STES.M.MX.PRMNTO01.IXNB.SA.2008100' %in% names(tts))
      stop("test series has disappeared from provider.")

if(! TSsdmx::verifyQuery('INEGI','DF_STES.M.MX.PRMNTO01.IXNB.SA.2008100'))
       stop("Query 2 does not verify. Provider changed something.")

tts[['DF_STES.M.MX.PRMNTO01.IXNB.SA.2008100']]

tts <- getSDMX("INEGI", 'DF_STES.M.MX.PRMNTO01.IXNB.SA.2008100') 
  
if ("Jan 1993" != start(tts[[1]]))  stop("INEGI test 2 start date error.")
if (12 != frequency(tts[[1]]))      stop("INEGI test 2 frequency error.")
if (1 != length(names(tts)))     stop("INEGI test 2 number of series changed.")


if (FALSE){ 

tts <- getSDMX("INEGI", 'DF_STES.M.MX.PRMNTO01.IXNB.SA.2008100', start='Feb 1995') 

# this part could use lots of cleanup, but it is slow and I am
# not sure if provider is still making changes.

#why slash? in DF_STEI/
  #tts = getTimeSeries("INEGI", "DF_STEI/..C1161+C1162+C5004.....") #slow, works sometimes
  tts = getSDMX("INEGI", "DF_STEI/..C1161+C1162+C5004.....") #slow, works sometimes
  
  nm <- getFlows('INEGI')  # can be slow
  length(nm)  # 8
  names(nm)
  getDimensions('INEGI','DF_COMTRADE') 
  getCodes('INEGI','DF_COMTRADE', 'FREQ') # fails sometimes
  names(getDimensions('INEGI','DF_COMTRADE')) 

  #if(! TSsdmx::verifyQuery('INEGI',
  #        'DF_COMTRADE.Q.MX.TOTAL.CAN.*.*.USD.Z.CAN.*'))
  #     stop("Query 1 does not verify. Provider changed something.")
  
 
  #####  FAILURE #####:  
 #13 BUG? or something else. Need working examples here 
  #tts <- getSDMX("INEGI", 'DF_COMTRADE.Q.MX.TOTAL.CAN.*.*.USD.Z.CAN.*') #empty
  tts <- getSDMX("INEGI", 'DF_COMTRADE.Q.MX.TOTAL.CAN.*.*.USD.*.*.*') #empty & slow
  names(tts)

#####  FAILURE #####:  
  tts <- getSDMX("INEGI", 'DF_COMTRADE.Q.MX.TOTAL.*.*.*.USD.*.*.*')  #empty & slow
  # (previously failed with: Comment must start with "<!--". )

####  FAILURE #####:  
  tts <- getSDMX("INEGI", 'DF_COMTRADE.Q.MX.TOTAL.*.*.*.*.*.*.*')  #empty
  tts <- getSDMX("INEGI", 'DF_COMTRADE.*.MX.TOTAL.*.*.*.*.*.*.*')  #empty
  
#####  FAILURE #####:  
  tts <- getSDMX("INEGI", 'DF_COMTRADE.*.*.TOTAL.*.*.*.*.*.*.*')  #empty 
  
#####  FAILURE #####:  
  tts <- getSDMX("INEGI", 'DF_COMTRADE.*.*.*.*.*.*.*.*.*.*') # very slow responding, then #Error in .jcall("RJavaTools", "Ljava/lang/Object;", "invokeMethod", cl,  : 
#    it.bankitalia.reri.sia.util.SdmxException: Exception. Class: java.net.SocketException #.Message: Connection reset
    
  #names(tts)
 
}
