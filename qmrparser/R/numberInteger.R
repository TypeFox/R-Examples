#do not edit, edit noweb/qmrparser.nw
numberInteger <- function(action = function(s) list(type="numberInteger",value=s), 
                          error = function(p) list(type="numberInteger",pos  =p))
  
  function (stream) {
    cstream <- streamParserNextChar(stream)
    
    if ( cstream$status == "eof" )   return(list(status="fail",node=error(streamParserPosition(stream)),stream=stream))
    
    if ( cstream$char %in% c("+","-") ) {
      s        <- cstream$char
      cstream  <- numberNatural()(cstream$stream) 
    }
    else {
      s       <- ""
      cstream <- numberNatural()(stream) 
    }   
    if ( cstream$status == "fail" )   return(list(status="fail",node=error(streamParserPosition(stream)),stream=stream))
    return(list(status="ok",node=action(paste(s,cstream$node$value,sep="")),stream=cstream$stream)) 
    
  }     

