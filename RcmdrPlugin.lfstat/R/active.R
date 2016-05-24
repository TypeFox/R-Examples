activelf <- function(){
  activeDataSetP() && inherits(justDoIt(activeDataSet()),'lfobj')
} 

  

activelfandbf <- function(){
  if( activelf()) {"baseflow" %in% Variables()}else(FALSE)
}

isthereanRFD <- function(){
  alle <- objects(envir=.GlobalEnv)
  b <- NULL
  for(ii in alle){
   b[ii] <- inherits(get(ii),"rfd")
  }
  any(b)
  }
    
