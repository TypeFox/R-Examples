#do not edit, edit noweb/qmrparser.nw
whitespace <- function(action = function(s) list(type="white",value=s),
                       error  = function(p) list(type="white",pos  =p) ) 
  
  function (stream) {
    
    cstream <- streamParserNextChar(stream)
    
    if ( cstream$status == "eof" || ! isWhitespace(cstream$char) ) return(list(status="ok",node=action(""),stream=stream)) 
    
    s <- cstream$char
    repeat {
      stream  <- cstream$stream
      cstream <- streamParserNextCharSeq(stream) 
      
      if (  cstream$status == "eof" || ! isWhitespace(cstream$char)  ) return(list(status="ok",node=action(paste(s,collapse="")),stream=stream)) 
      
      s      <- c(s,cstream$char)
    }
  }

