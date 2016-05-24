#do not edit, edit noweb/qmrparser.nw
eofMark  <- function(action = function(s) list(type="eofMark",value=s),
                     error  = function(p) list(type="eofMark",pos=p  )) 

  function (stream) {
                  
    cstream <- streamParserNextChar(stream)
    
    if ( cstream$status != "eof" ) return(list(status="fail",node=error(streamParserPosition(stream)), stream=stream))
    
    return(list(status="ok",node=action(""),stream=cstream$stream))
  }
