#do not edit, edit noweb/qmrparser.nw
concatenation <- function(..., 
                          action = function(s)  list(type="concatenation", value=s),
                          error  = function(p,h) list(type="concatenation", pos=p   ,h=h)) 
  
  function(stream) {
    
    streamFail <- stream
    dots       <- list(...)
    value      <- rep(list(list()),length(dots))
    
    for( i in 1:length(dots) ) {        
      cstream     <- (dots[[i]])(stream)
      
      if (cstream$status == "fail") return(list(status="fail",node=error(streamParserPosition(streamFail),cstream$node),stream=streamFail))
      
      value[[i]]  <- cstream$node      
      stream      <- cstream$stream
    }
    return(list(status="ok",node=action(value),stream=stream))
##
  }

