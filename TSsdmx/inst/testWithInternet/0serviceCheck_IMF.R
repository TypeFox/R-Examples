######################## IMF #######################
require("RJSDMX")

#  sdmxHelp()
  
  if(! 'PGI' %in% names(getFlows('IMF')))
     stop("PGI has disappearde from IMF flows (again). Provider changed something.")

  # IMF PGI codes were not working for some period prior to Feb 9, 2015. 
  #  but then started working again. 

  # PGI
  #  REF_AREA.DATASOURCE.PGI_CONCEPT.FREQ.UNIT0FMEASURE
  #     CA: Canada
  #  INDICATOR
  #    003: National Accounts
  #    IFS: International Financial Statistics
  #    AIP: Industrial Production
  #    BIS_BP6: Balance on Secondary Income
  #    NCG: Government Consumption Expenditure
  #    NGDP: Gross Domestic Product (Nominal)
  #  DATA SOURCE
  #    PGI: Principal Global Indicators
  #  UNIT
  #    L:  USD  
  #    N:   National currency 
  #    NSA: National currency SA
  #  FREQ
  #    A:  
  #    M:  
  #    Q:  

  names(getDimensions('IMF','PGI')) 
  getCodes('IMF','PGI', 'FREQ')
  
  if(! TSsdmx::verifyQuery('IMF', 'PGI.CA.*.*.*.*'))
     stop("Query 1 does not verify. Provider changed something.")
  
  tts0 <- getSDMX('IMF', 'PGI.CA.*.*.*.*')   # length # 627   #774 Feb 9, 2015
  nm <- names(tts0)
  length(nm) 
  
  nm[grepl('PGI.CA.BIS.', nm )] 
  #[1] "PGI.CA.BIS_BP6.PGI.L.A" "PGI.CA.BIS_BP6.PGI.L.Q"  # Feb 9, 2015

  # note that grepl uses . as any char so this gets above
  #z <- tts0[grepl('PGI.CA.BIS.', nm )]
    
  z <- getSDMX('IMF', 'PGI.CA.BIS_BP6.PGI.L.A')
  # start was 2005 for awhile (circa spring 2015)
  if(start(z[[1]]) !=  1948)  stop("test 1 start date changed (again).")
  if(frequency(z[[1]]) !=  1) stop("test 1  frequency changed.")

  tts <- getSDMX('IMF', 'PGI.CA.BIS_BP6.*.L.Q')	
  names(tts)
  
  #  at one time this was TRUE, but not Feb 9, 2015 
  #	"PGI.CA.BIS.FOSAB.Q.L_M" %in% nm 

  if( TSsdmx::verifyQuery('IMF', 'PGI.CA.BIS.*.*.*'))
     stop("Query 2 now verifies. Provider added it again.")

  if( TSsdmx::verifyQuery('IMF', 'PGI.CA.IFS.*.*.*'))
     stop("Query 3 now verifies. Provider added it again.")
 
  #tts <- getSDMX('IMF', 'PGI.CA.IFS.*.Q.N_M') #fails (empty result)

