if (FALSE){ 
######################## National Bank Belguim #######################

  require("TSsdmx")
  
  nbb <- TSconnect("sdmx", dbname="NBB")

##  BUG
## if( ! verifyQuery('NBB', 'HICP.000000.*.*', verbose=FALSE))
##     stop("verifyQuery NBB wildcards failed")

######### annual

  z <- TSget('HICP.000000.BE.A', nbb) 

  if(!all(start(z) ==  c(1992,1))) stop("test 1 start date changed.")
  if(frequency(z) !=  1)           stop("test 1  frequency changed.")
  if(tframe::nseries(z) !=  1)     stop("test 1  number of series changed.")

  z <- TSget('HICP.*.BE.A', start=c(2000,1), end=c(2010,1), nbb) 

  if(!all(start(z) ==  c(2000,1)))  stop("test 2 start date changed.")
  if(!all(end(z)   ==  c(2010,1)))  stop("test 2  end  date changed.")
  if(frequency(z) !=  1)            stop("test 2  frequency changed.")
  if(tframe::nseries(z) !=  9)      stop("test 2  number of series changed.")

######### monthly

  z <- TSget('HICP.000000.BE.M', nbb) 

  if(!all(start(z) ==  c(1992,1)))  stop("test 3 start date changed.")
  if(frequency(z) !=  12)           stop("test 3  frequency changed.")
  if(tframe::nseries(z) !=  1)      stop("test 3  number of series changed.")

  z <- TSget('HICP.*.BE.M', start=c(2000,1), end=c(2010,1), nbb) 

  if(!all(start(z) ==  c(2000,1)))  stop("test 4 start date changed.")
  if(!all(end(z)   ==  c(2010,1)))  stop("test 4  end  date changed.")
  if(frequency(z) !=  12)           stop("test 4  frequency changed.")
  if(tframe::nseries(z) !=  9)      stop("test 4  number of series changed.")

}
