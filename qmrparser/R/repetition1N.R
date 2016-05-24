#do not edit, edit noweb/qmrparser.nw
repetition1N <- function(rpa, 
                         action = function(s)   list(type="repetition1N",value=s  ),
                         error  = function(p,h) list(type="repetition1N",pos=p,h=h))

  
  function(stream) {
    cstream <- rpa(stream)
    
    if ( cstream$status == "fail" )  return(list(status="fail",node=error(streamParserPosition(stream),cstream$node),stream=stream))
    
    iinc        <- 1000
    value       <- rep(list(list()),iinc)
    imax        <- iinc
    i           <- 1    
    value[[i]]  <- cstream$node    
    stream      <- cstream$stream 
    
    while( TRUE ) {
      cstream <- rpa(stream)
      if ( cstream$status == "fail" ) break()
      i           <- i + 1
      if ( i >= imax ) { 
        value <- c(value, rep(list(list()),iinc))
        imax  <- imax + iinc
      }        
      value[[i]]  <- cstream$node
      stream      <- cstream$stream
    }
    return(list(status="ok",node=action(value[1:i]),stream=stream))
  }
