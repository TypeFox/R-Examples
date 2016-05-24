#do not edit, edit noweb/qmrparser.nw
string  <- function(isQuote= function(c) switch(c,'"'=,"'"=TRUE,FALSE),
                    action = function(s) list(type="string",value=s),
                    error  = function(p) list(type="string",pos  =p)) 
  
  function (stream) {
    cstream   <- streamParserNextChar(stream)
    
    if ( cstream$status == "eof" || ! isQuote(cstream$char) )      return(list(status="fail",node=error(streamParserPosition(stream)),stream=stream))  
    
    delimiter <- cstream$char
    
    charPrev <- ""
    s        <- ""
    repeat {
      cstream <- streamParserNextCharSeq(cstream$stream)
      
      if ( cstream$status == "eof" )   return(list(status="fail",node=error(streamParserPosition(stream)),stream=stream))  
      
      if ( cstream$char == delimiter && charPrev != "\\" )  return(list(status="ok"  ,node=action(paste(s,collapse="")), stream=cstream$stream))
      
      charPrev <- cstream$char
      s        <- c(s,cstream$char)
    }
  } 
