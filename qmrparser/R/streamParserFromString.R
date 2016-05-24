#do not edit, edit noweb/qmrparser.nw
streamParserFromString <- function(string) {
  
  stringIn <- paste(string,collapse="\n")
  nchars   <- nchar(stringIn)
  
  ## One at a time character processing multiplies space by 8 (in 64 bits OS) but reduces time /40
  stringIn <- strsplit(stringIn,split=character(0))[[1]]
  string   <- stringIn
  if( nchars > 0 && length(string) < 1 ) stop("Error: encoding error")
  
  return(
         list(
              streamParserNextChar =  function(stream) {
                if ( stream$pos + 1 > stream$lenchar ) list(status="eof",char=""  ,stream=stream)
                else { 
                  
                  stream$pos     <- stream$pos+1
                  
## version for substr char           <- substr(string, stream$pos, stream$pos)                  
                  char           <- string[stream$pos]                  
                  
                  if ( char == "\n" ) { 
                    stream$line    <- stream$line + 1
                    stream$linePos <- 0              
                  } else            { 
                    stream$linePos <- stream$linePos + 1
                  }
                  
                  list(status="ok" ,char=char,stream=stream)
                }
              },
              
              streamParserNextCharSeq =  function(stream) {
                if ( stream$pos + 1 > stream$lenchar ) list(status="eof",char=""  ,stream=stream)
                else { 
                  
                  stream$pos     <- stream$pos+1

## version for substr char <- substr(string,stream$pos, stream$pos)                                    
                  char           <- string[stream$pos]                  
                  
                  if ( char == "\n" ) { 
                    stream$line    <- stream$line + 1
                    stream$linePos <- 0              
                  } else            { 
                    stream$linePos <- stream$linePos + 1
                  }
                  
                  list(status="ok" ,char=char,stream=stream)
                }
              },
              
              streamParserClose      =  function(stream) { lenchar <- 0 ; invisible(NULL) 
                                                        },
              streamParserPosition  =   function(stream) { list(fileName="", line=stream$line, linePos=stream$linePos+1, streamPos=stream$pos+1 )
                                                        },
              
##            lenchar  = nchar(stringIn) , ## version for substr
              lenchar  = length(stringIn),
              pos      = 0,
              line     = 1,
              linePos  = 0
              )
         )
}
