#do not edit, edit noweb/qmrparser.nw
commentParser <- function(beginComment,endComment,
                    action = function(s)list(type="commentParser",value=s),
                    error  = function(p)list(type="commentParser",pos=p  ))
  
  function (stream) {

    ## begin comment    
    s       <- as.character()
    cstream <- streamParserNextChar(stream)             
    
    for( i in 1:nchar(beginComment) ) {
      
      if ( cstream$status == "eof" || cstream$char != substr(beginComment,i,i) ) return(list(status="fail",node=error(streamParserPosition(stream)),stream=stream))

      s       <- c(s,cstream$char)
      
      cstream <- streamParserNextCharSeq(cstream$stream)       
    }

    previo        <- ""
    lenEndComment <- nchar(endComment)
    repeat {    
      ## end comment     
      init.s       <- s
      init.cstream <- cstream      
      for( i in 1:lenEndComment ) {
        
        if ( cstream$status == "eof" ) return(list(status="fail",node=error(streamParserPosition(stream)),stream=stream))

        s       <- c(s,cstream$char)
        
        if ( cstream$char != substr(endComment,i,i) ) break() 
        if ( previo       == "\\"                   ) break() 
      
        if ( i == lenEndComment) return(list(status="ok", node=action(paste(s,collapse="")), stream=cstream$stream)) 
        
        cstream <- streamParserNextCharSeq(cstream$stream)           
      }
        previo  <- init.cstream$char
        s       <- c(init.s,previo)        
        cstream <- streamParserNextCharSeq(init.cstream$stream)         
      
    }

  }

