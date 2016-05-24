#do not edit, edit noweb/qmrparser.nw
symbolic <- function(charFirst=isLetter,
                     charRest=function(ch) isLetter(ch) || isDigit(ch) || ch == "-",
                     action = function(s) list(type="symbolic",value=s), 
                     error  = function(p) list(type="symbolic",pos  =p)) 
  
  function (stream) {
    cstream <- streamParserNextChar(stream)             
    
    if ( cstream$status == "eof" || ! charFirst(cstream$char) ) return(list(status="fail",node=error(streamParserPosition(stream)),stream=stream))     
    
    s <- cstream$char
    repeat {
      stream  <- cstream$stream
      cstream <- streamParserNextCharSeq(stream) 
      
      if (  cstream$status == "eof" || ! charRest(cstream$char) )  return(list(status="ok",node=action(paste(s,collapse="")),stream=stream)) 

      s      <- c(s,cstream$char)
    }
  }

