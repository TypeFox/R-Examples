#do not edit, edit noweb/qmrparser.nw
numberScientific <- function(action = function(s) list(type="numberScientific",value=s), 
                             error  = function(p) list(type="numberScientific",pos=p  )
                             )
  function (stream) {
    streamFail <- stream
    
####                    numberFloat
    cstream    <- streamParserNextChar(stream)
    
    if ( cstream$status == "eof" )   return(list(status="fail",node=error(streamParserPosition(streamFail)),stream=streamFail))
    
    if( cstream$char == "-" || cstream$char == "+" ){ 
      signo   <- cstream$char 
      stream  <- cstream$stream
      cstream <- streamParserNextChar(stream) 
      
      if ( cstream$status == "eof" )   return(list(status="fail",node=error(streamParserPosition(streamFail)),stream=streamFail))
    } else {
      signo <- ""
    }
    if ( cstream$char == '.' )  {
      entero   <- ""
      punto    <- cstream$char
      
#####                           numberNatural
      cstream <- streamParserNextChar(cstream$stream)
      
      if ( cstream$status == "eof" || ! isDigit(cstream$char) ) return(list(status="fail",node=error(streamParserPosition(streamFail)),stream=streamFail)) 
      
      s <- cstream$char
      repeat {
        stream  <- cstream$stream
        cstream <- streamParserNextCharSeq(stream) 
        if ( cstream$status == "eof" || ! isDigit(cstream$char) ) break()
        s       <- c(s,cstream$char)
      }
##### 
      decimal <- paste(s,collapse="")
    }
    else {
#####                           numberNatural
      if ( ! isDigit(cstream$char) ) return(list(status="fail",node=error(streamParserPosition(streamFail)),stream=streamFail)) 
      s <- cstream$char
      repeat {
        stream  <- cstream$stream
        cstream <- streamParserNextCharSeq(stream) 
        if ( cstream$status == "eof" || ! isDigit(cstream$char) ) break()
        s       <- c(s,cstream$char)
      }
#####
      entero  <- paste(s,collapse="")
##
      cstream <- streamParserNextChar(stream)
      if ( cstream$char == '.' )  {
        punto    <- cstream$char
        stream   <- cstream$stream
#####                           numberNatural
        cstream <- streamParserNextChar(stream)
        if ( cstream$status == "eof" || ! isDigit(cstream$char) ) {
          decimal <- "" 
        } else {
          s <- cstream$char
          repeat {
            stream  <- cstream$stream
            cstream <- streamParserNextCharSeq(stream) 
            if ( cstream$status == "eof" || ! isDigit(cstream$char) ) break()
            s      <- c(s,cstream$char)
          }
          decimal  <- paste(s,collapse="")
        }
#####
      } else {
        punto   <- ""
        decimal <- ""
      } 
    }
    mantisa <- paste(signo,entero,punto,decimal,sep="")
#
    cstream <- streamParserNextChar(stream) 
    if ( cstream$char == "E" || cstream$char == "e" ) { 
      E        <- cstream$char
      stream   <- cstream$stream
      cstream  <- streamParserNextChar(stream)
      
      if ( cstream$status == "eof" ) return(list(status="fail",node=error(streamParserPosition(streamFail)),stream=streamFail))
#####
      if ( cstream$char == "-" || cstream$char == "+" ) { 
        signoE  <- cstream$char 
        stream  <- cstream$stream
        cstream <- streamParserNextChar(stream) 
        
        if ( cstream$status == "eof" )   return(list(status="fail",node=error(streamParserPosition(streamFail)),stream=streamFail))
      } else {
        signoE <- ""
      }
#####                           numberNatural
      if ( ! isDigit(cstream$char) ) return(list(status="fail",node=error(streamParserPosition(streamFail)),stream=streamFail)) 
      s <- cstream$char
      repeat {
        stream  <- cstream$stream
        cstream <- streamParserNextCharSeq(stream) 
        if ( cstream$status == "eof" || ! isDigit(cstream$char) ) break()
        s       <- c(s,cstream$char)
      }
#####
      exponente  <- paste(signoE,paste(s,collapse=""),sep="")
    } else {
      E          <- ""
      exponente  <- ""
      
    }
    return(list(status="ok",node=action(paste(mantisa,E,exponente,sep="")),stream=stream))
  }
