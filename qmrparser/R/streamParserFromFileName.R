#do not edit, edit noweb/qmrparser.nw
streamParserFromFileName <- function(fileName,encoding = getOption("encoding")) {

  ## streamParseFromFileName can be used?
  if( Sys.info()["sysname"] == "Windows" ) {
    
    fromString <- TRUE 
    
  } else {
    conn <- file(fileName,"r",encoding =encoding)  
    if ( ! isOpen(conn)     ) stop(paste("Error: file cannot be opened",fileName)) 
    fromString <- tryCatch({ seek(conn) ; FALSE}, error =function(e) TRUE, finally= close(conn) )
    
  }
  
  if ( fromString ) return( streamParserFromString( readLines( fileName, encoding=encoding)) )
  else
    return( list(
                 streamParserNextChar =  function(stream) {
                   
                   if ( stream$pos != seek(stream$conn) ) seek(stream$conn,stream$pos)                 
                                      
                   char <- readChar(stream$conn,nchars=1,useBytes = FALSE)
                   
                   if (length(char) == 0)       
                     list(status="eof",char=""  ,stream=stream) 
                   else {
                     
                     stream$pos    <- seek(stream$conn)           
                                          
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
                   char <- readChar(stream$conn,nchars=1,useBytes = FALSE)
                   if (length(char) == 0)       
                     list(status="eof",char=""  ,stream=stream)
                   
                   else {
                     
                     stream$pos       <- seek(stream$conn)                                        
                     
                     if ( char == "\n" ) { 
                       stream$line    <- stream$line + 1
                       stream$linePos <- 0              
                     } else            { 
                       stream$linePos <- stream$linePos + 1
                     }
                     
                     list(status="ok" ,char=char,stream=stream)
                   }
                 },
                 
                 streamParserClose      =  function(stream) { close(stream$conn) ; stream$conn <- -1 ; invisible(NULL) 
                                                            },
                 
                 streamParserPosition  =   function(stream) { list(fileName=stream$fileName, line=stream$line, linePos=stream$linePos+1, streamPos=stream$pos+1)
                                                            },
                 conn     = local({ 
                   conn <- file(fileName,"r",encoding =encoding)  
                   if ( ! isOpen(conn)     ) stop(paste("Error: file cannot be opened.",fileName)) 
                   tryCatch( seek(conn) , error =function(e) stop(paste("Error:  'seek' not enabled for this connection", fileName)))                                        
                   conn 
                 }),
                 pos      = 0,
                 line     = 1,
                 linePos  = 0,
                 fileName = fileName
                 )
           )
}
