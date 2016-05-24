if (FALSE){ # previously this all worked
########################### "INEGI" ###########################
##### Instituto Nacional de Estadistica y Geografia (Mexico) ######

require("TSsdmx")
  
inegi <- TSconnect("sdmx", dbname="INEGI")

z <- TSget('DF_STES.M.MX.*.*.*.*', inegi) 
  
if(!all(start(z) ==  c(1991,1))) stop("test 1 start date changed.")
if (12 != frequency(z) )         stop("test 1 frequency error.")
if(tframe::nseries(z) !=  26)    stop("test 1  number of series changed.")

z <- TSget('DF_STES.M.MX.PRMNTO01.IXNB.SA.2008100', inegi) 
  
if(!all(start(z) ==  c(1993,1))) stop("test 2 start date changed.")
if (12 != frequency(z))          stop("test 2 frequency error.")
if(tframe::nseries(z) !=  1)     stop("test 2  number of series changed.")

}
