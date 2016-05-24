#do not edit, edit noweb/qmrparser.nw
separator <- function(action= function(s) list(type="separator",value=s), 
                      error = function(p) list(type="separator",pos=p  )) 
  function (stream) {
    cstream <- streamParserNextChar(stream)
    
    if ( cstream$status == "eof" ) return(list(status="fail",node=error(streamParserPosition(stream)),stream=stream))
    
    ## possible whites
    s <- ""
    repeat {
      if (  cstream$status == "eof" || ! isWhitespace(cstream$char)  ) break() 
      s       <- c(s,cstream$char)
      stream  <- cstream$stream
      cstream <- streamParserNextCharSeq(stream) 
    }
    
    ##possible , ;
    if ( cstream$status == "eof" || ( cstream$char != ',' && cstream$char != ';') ) {
      if ( length(s) > 1 ) return(list(status="ok",node=action(paste(s,collapse="")),stream=stream))
      else
        return(list(status="fail",node=error(streamParserPosition(stream)),stream=stream))
    }
    ##    
    s       <- c(s,cstream$char)
    
    ## possible white
    repeat {
      stream  <- cstream$stream
      cstream <- streamParserNextCharSeq(stream) 
      
      if (  cstream$status == "eof" || ! isWhitespace(cstream$char)  ) return(list(status="ok",node=action(paste(s,collapse="")),stream=stream)) 
      
      s       <- c(s,cstream$char)
    }
  }
