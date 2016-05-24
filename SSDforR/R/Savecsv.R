Savecsv<-
  function(){
    
  
    outdat <- get("ssd", envir  = environment())
    if (exists('ssd')) write.csv(outdat,file = tclvalue(tcl("tk_getSaveFile")),row.names=FALSE)
  }              