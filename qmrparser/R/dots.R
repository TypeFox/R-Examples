#do not edit, edit noweb/qmrparser.nw
dots <- function(action = function(s) list(type="dots",value=s), 
                 error  = function(p) list(type="dots",pos  =p)) 
  
  function (stream) {
    cstream <- streamParserNextChar(stream)             
    
    if ( cstream$status == "eof" || cstream$char != "." ) return(list(status="fail",node=error(streamParserPosition(stream)),stream=stream)) 
    s <- cstream$char
    repeat {
      stream  <- cstream$stream
      cstream <- streamParserNextCharSeq(stream) 
      
      if ( cstream$status == "eof" || ! ( cstream$char == "." ) ) return(list(status="ok",node=action(paste(s,collapse="")),stream=stream)) 
      
      s      <- c(s,cstream$char)
    }
  }
