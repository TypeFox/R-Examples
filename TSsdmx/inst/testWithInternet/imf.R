#http://www.imf.org

require("TSsdmx")

# RJSDMX::sdmxHelp()  # can be useful for finding series identifiers, etc

imf <- TSconnect("sdmx",  dbname="IMF")

 if(FALSE != verifyQuery('IMFx', 'PGI.CA.*.*.*.*', verbose=FALSE)) 
     stop("verifyQuery provider IMFx failed")

 if( ! verifyQuery('IMF',  'PGI.CA.*.*.*.*', verbose=FALSE))
     stop("verifyQuery IMF wildcards failed")

 if(FALSE != verifyQuery('IMF', 'PGI.CAN.*.*.*.*', verbose=FALSE)) 
     stop("verifyQuery bad dimension check failed")
 
if (FALSE){
  
  z <- TSget("PGI.CA.BIS.FOSLB.A.L_M", imf)  #13 BUG? this was previously not empty
   
#####  FAILURE #####:   empty result but retrieved above
  z <- TSget('PGI.CA.BIS.*.*.L_M', imf) #fails (empty result)
  
#####  FAILURE #####:   empty result but retrieved above
  z <- TSget("PGI.CA.BIS.FOSAB.Q.L_M") #fails (empty result)

  #  even though it was returned above
  	"PGI.CA.BIS.FOSAB.Q.L_M" %in% seriesNames(z)
  #   and 

  nm <- seriesNames(z)
  
  nm["PGI.CA.BIS.FOSAB.Q.L_M" ]

 
  nm[grepl('PGI.CA.IFS.', nm )] # this suggests these should work but

#####  FAILURE #####:   empty result but retrieved above
  tts <- getSDMX('IMF', 'PGI.CA.IFS.*.Q.N_M') #fails (empty result)
  names(tts)

}
