#do not edit, edit noweb/qmrparser.nw
alternation <- function(..., 
                        action = function(s)     list(type="alternation",value=s), 
                        error  = function(p,h) list(type="alternation",pos  =p,h=h) ) 
  
  function(stream) {
    dots       <- list(...)
    h          <- rep(list(list()),length(dots))
    
    for( i in 1:length(dots) ) {        
      cstream <- (dots[[i]])(stream)
      
      if (cstream$status=="ok") return(list(status="ok",node=action(cstream$node),stream=cstream$stream))
      
      h[[i]] <- cstream$node
    }
    return(list(status="fail",node=error(streamParserPosition(stream),h),stream=stream))
    
  }
