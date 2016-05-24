#do not edit, edit noweb/qmrparser.nw
charInSetParser <- function(fun, 
                           action = function(s) list(type="charInSet",value=s), 
                           error  = function(p) list(type="charInSet",pos  =p)) 
  
  function(stream) {
    cstream <- streamParserNextChar(stream)
    
    if ( cstream$status == "eof" )   return(list(status="fail",node=error(streamParserPosition(stream)),stream=stream))

    if ( fun(cstream$char) )    
      return(list(status="ok", node=action(cstream$char), stream=cstream$stream))
    return(list(status="fail",node=error(streamParserPosition(stream)),stream=stream))                                 
  }
