#do not edit, edit noweb/qmrparser.nw
keyword <- function(word,
                    action = function(s) list(type="keyword",value=s), 
                    error  = function(p) list(type="keyword",pos  =p)) 

  function (stream) {

    ## begin word
    s       <- as.character()
    cstream <- streamParserNextChar(stream)             
    
    for( i in 1:nchar(word) ) {
          
      if ( cstream$status == "eof" || cstream$char != substr(word,i,i) ) return(list(status="fail",node=error(streamParserPosition(stream)),stream=stream))

      s       <- c(s,cstream$char)

      if( i == nchar(word) )
      return(list(status="ok", node=action(paste(s,collapse="")), stream=cstream$stream)) 
      
      cstream <- streamParserNextChar(cstream$stream)                 
      
    }
    
  }
