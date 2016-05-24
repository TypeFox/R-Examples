#do not edit, edit noweb/qmrparser.nw
numberFloat <- function(action = function(s) list(type="numberFloat",value=s), 
                        error  = function(p) list(type="numberFloat",pos=p) ) 
  
  function (stream) {
    streamFail <- stream
    cstream    <- streamParserNextChar(stream)
    
    if ( cstream$status == "eof" )   return(list(status="fail",node=error(streamParserPosition(streamFail)),stream=streamFail))

    ## sign
    if ( cstream$char %in% c("+","-") ) { 
      signo   <- cstream$char 
      stream  <- cstream$stream
      cstream <- streamParserNextChar(stream) 
      
      if ( cstream$status == "eof" )   return(list(status="fail",node=error(streamParserPosition(streamFail)),stream=streamFail))
      
    } else {
      signo <- ""
    }
    if ( cstream$char == '.' )  {
      entero   <- ""
      punto    <- cstream$char
      cstream  <- numberNatural()(cstream$stream) 
      if ( cstream$status == "ok" ) {
        decimal <- cstream$node$value
        stream  <- cstream$stream
      } else {
        return(list(status="fail",node=error(streamParserPosition(streamFail)),stream=streamFail))
        
      }                 
    }
    else {
      cstream <- numberNatural()(stream)
      
      if (cstream$status=="fail") return(list(status="fail",node=error(streamParserPosition(streamFail)),stream=streamFail))
      
      entero  <- cstream$node$value
      stream  <- cstream$stream
      cstream <- streamParserNextChar(stream)
      if ( cstream$char == '.' )  {
        punto    <- cstream$char
        stream   <- cstream$stream
        cstream  <- numberNatural()(stream)
        if ( cstream$status == "ok" ) {
          decimal <- cstream$node$value
          stream  <- cstream$stream
        } else {
          decimal <- ""
        }
      } else {
        punto   <- ""
        decimal <- ""
      } 
    }
    return(list(status="ok",node=action(paste(signo,entero,punto,decimal,sep="")),stream=stream)) 
    
  }
