#do not edit, edit noweb/qmrparser.nw
numberNatural <- function(action = function(s) list(type="numberNatural",value=s), 
                          error  = function(p) list(type="numberNatural",pos  =p) ) 
  
  function (stream) {
    cstream <- streamParserNextChar(stream)
    
    if ( cstream$status == "eof" || ! isDigit(cstream$char) ) return(list(status="fail",node=error(streamParserPosition(stream)),stream=stream))
    
    s <- cstream$char
    
    repeat {
      
      stream  <- cstream$stream
      cstream <- streamParserNextCharSeq(stream) 
      
      if ( cstream$status == "eof" || ! isDigit(cstream$char) )          return(list(status="ok",node=action(paste(s,collapse="")),stream=stream)) 
      
      s      <- c(s,cstream$char)
      
    }
  }
