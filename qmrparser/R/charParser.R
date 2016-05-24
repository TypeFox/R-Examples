#do not edit, edit noweb/qmrparser.nw
charParser <- function(char, 
                      action = function(s) list(type="char",value=s), 
                      error  = function(p) list(type="char",pos  = p)) 
  
  function(stream) {
    cstream <- streamParserNextChar(stream)
    
    if ( cstream$status == "eof" )   return(list(status="fail",node=error(streamParserPosition(stream)),stream=stream))
    
    if ( cstream$char   != char  )   return(list(status="fail",node=error(streamParserPosition(stream)),stream=stream))  
    return(list(status="ok"  ,node=action(char),stream=cstream$stream))
  }
##
