  printCStream <- function(cstream) print(paste("[",cstream$node$value,"] [",cstream$status,"]",sep=""))
  checkTokenParserOk  <- function(text,parse,type,token=text,msg=paste(type,": [",token,"]",sep=""))
                  checkIdentical(list(status="ok",node=list(type=type ,value=token)),parse(streamParserFromString(text)) [c("status","node")],msg=msg)
  checkTokenParserFail <- function(text,parse,type,token=text,msg=paste(type,": [",token,"]",sep=""))
                          checkIdentical("fail",parse(streamParserFromString(text))$status,msg=msg)
  checkSentParserOk <- function(text,parse,type,sent=text,msg=paste(type,": [",sent,"]",sep="")) {
                                  r            <- parse(streamParserFromString(text)) [c("status","node")]
                                  r$node$value <- paste(r$node$value,collapse="",sep="")
  #print(r)
                                  checkIdentical(list(status="ok",node=list(type=type ,value=sent)),r,msg=msg)
  }
  checkSentParserFail  <- function(text,parse,type,sent=text,msg=paste(type,": [",sent,"]",sep="")) {
                  checkIdentical("fail",parse(streamParserFromString(text))$status,msg=msg)
  }
  action       <- function(s) paste(s,collapse="",sep="")
  actionString <- function(s) paste("'",s,"'",collapse="",sep="")
test.tokens01 <- function() {
  checkTokenParserOk     ("",eofMark(),"eofMark")
  checkTokenParserOk  (""  ,whitespace(),"white")
  checkTokenParserOk  ("  ",whitespace(),"white")
  checkTokenParserOk  ("\n",whitespace(),"white")
  checkTokenParserOk  ("A" ,whitespace(),"white","")
  checkTokenParserOk  ("1" ,whitespace(),"white","")
  checkTokenParserFail (""       ,separator(),"separator")
  checkTokenParserFail ("A"      ,separator(),"separator")
  checkTokenParserFail ("1"      ,separator(),"separator")
  checkTokenParserOk   (" "      ,separator(),"separator")
  checkTokenParserOk   ("\n"     ,separator(),"separator")
  checkTokenParserOk   (","      ,separator(),"separator")
  checkTokenParserOk   (";"      ,separator(),"separator")
  checkTokenParserOk   (" ;"     ,separator(),"separator")
  checkTokenParserOk   ("; "     ,separator(),"separator")
  checkTokenParserOk   (" ; "    ,separator(),"separator")
  checkTokenParserOk   ("  ,  "  ,separator(),"separator")
  checkTokenParserOk  ("123"       ,numberNatural()   ,"numberNatural"    )
  checkTokenParserOk  ("123"       ,numberInteger()   ,"numberInteger"    )
  checkTokenParserOk  ("+123"      ,numberInteger()   ,"numberInteger"    )
  checkTokenParserOk  ("-123"      ,numberInteger()   ,"numberInteger"    )
  checkTokenParserOk  ("123"       ,numberFloat()     ,"numberFloat"      )
  checkTokenParserOk  ("+123"      ,numberFloat()     ,"numberFloat"      )
  checkTokenParserOk  ("-123"      ,numberFloat()     ,"numberFloat"      )
  checkTokenParserOk  ("123."      ,numberFloat()     ,"numberFloat"      )
  checkTokenParserOk  ("+123."     ,numberFloat()     ,"numberFloat"      )
  checkTokenParserOk  ("-123."     ,numberFloat()     ,"numberFloat"      )
  checkTokenParserOk  (".123"      ,numberFloat()     ,"numberFloat"      )
  checkTokenParserOk  ("+.123"     ,numberFloat()     ,"numberFloat"      )
  checkTokenParserOk  ("-.123"     ,numberFloat()     ,"numberFloat"      )
  checkTokenParserOk  ("127"          ,numberScientific(),"numberScientific" )
  checkTokenParserOk  ("+127"         ,numberScientific(),"numberScientific" )
  checkTokenParserOk  ("-127"         ,numberScientific(),"numberScientific" )

  checkTokenParserOk  ("127."          ,numberScientific(),"numberScientific" )
  checkTokenParserOk  ("+127."         ,numberScientific(),"numberScientific" )
  checkTokenParserOk  ("-127."         ,numberScientific(),"numberScientific" )

  checkTokenParserOk  (".127"          ,numberScientific(),"numberScientific" )
  checkTokenParserOk  ("+.127"         ,numberScientific(),"numberScientific" )
  checkTokenParserOk  ("-.127"         ,numberScientific(),"numberScientific" )

  checkTokenParserOk  ("0123.123"     ,numberScientific(),"numberScientific" )
  checkTokenParserOk  ("+0123.123"     ,numberScientific(),"numberScientific" )
  checkTokenParserOk  ("-0123.123"     ,numberScientific(),"numberScientific" )

  checkTokenParserOk  ("127E33"          ,numberScientific(),"numberScientific" )
  checkTokenParserOk  ("+127E33"         ,numberScientific(),"numberScientific" )
  checkTokenParserOk  ("-127E33"         ,numberScientific(),"numberScientific" )

  checkTokenParserOk  ("127.E33"          ,numberScientific(),"numberScientific" )
  checkTokenParserOk  ("+127.E33"         ,numberScientific(),"numberScientific" )
  checkTokenParserOk  ("-127.E33"         ,numberScientific(),"numberScientific" )

  checkTokenParserOk  (".127E33"          ,numberScientific(),"numberScientific" )
  checkTokenParserOk  ("+.127E33"         ,numberScientific(),"numberScientific" )
  checkTokenParserOk  ("-.127E33"         ,numberScientific(),"numberScientific" )

  checkTokenParserOk  ("0123.123E33"     ,numberScientific(),"numberScientific" )
  checkTokenParserOk  ("+0123.123E33"     ,numberScientific(),"numberScientific" )
  checkTokenParserOk  ("-0123.123E33"     ,numberScientific(),"numberScientific" )

  checkTokenParserOk  ("-121.01e-222" ,numberScientific(),"numberScientific" )
  checkTokenParserOk  ("1211e+33"      ,numberScientific(),"numberScientific" )
  checkTokenParserOk  ("-123e-222"    ,numberScientific(),"numberScientific" )

  checkTokenParserOk  ("+0123.123e+11"      ,numberScientific(),"numberScientific" )
  checkTokenParserOk  ("-1233.e-66"         ,numberScientific(),"numberScientific" )
  checkTokenParserOk  ("-1233.e+66"         ,numberScientific(),"numberScientific" )
  checkTokenParserOk  ("+.123E600"          ,numberScientific(),"numberScientific" )

  checkTokenParserFail("+.E600"             ,numberScientific(),"numberScientific" )
  checkTokenParserFail(".E600"              ,numberScientific(),"numberScientific" )
  checkTokenParserOk  ("abd"       ,symbolic(),"symbolic")
  checkTokenParserOk  ("ac12"      ,symbolic(),"symbolic")
  checkTokenParserOk  ("asd:"      ,symbolic(),"symbolic","asd")
  checkTokenParserOk  ('"aaaa"'    ,string(),"string","aaaa")
  checkTokenParserOk  ("'aaaa'"    ,string(),"string","aaaa")
  checkTokenParserOk  ('"aaa4567"' ,string(),"string","aaa4567")
  checkTokenParserOk  ("'aa aa'"   ,string(),"string","aa aa")
  checkTokenParserOk  ("'aa\\' aa'",string(),"string","aa\\' aa")
  checkTokenParserOk  ("?aaaa?"    ,string(isQuote=function(c) c=="?"),"string","aaaa")
    
  checkTokenParserOk  ("(**)"         ,commentParser("(*","*)"), "commentParser")
  checkTokenParserOk  ("(* *)"        ,commentParser("(*","*)"), "commentParser")
  checkTokenParserOk  ("(* 123 *)"    ,commentParser("(*","*)"), "commentParser")
  checkTokenParserOk  ("(* (*  *)"    ,commentParser("(*","*)"), "commentParser")
  checkTokenParserOk  ("(* (* * ) *)" ,commentParser("(*","*)"), "commentParser")
  checkTokenParserOk  ("/* (* * ) */" ,commentParser("/*","*/"), "commentParser")
  checkTokenParserOk  ("-- aadad\n"   ,commentParser("--","\n"), "commentParser")

  checkTokenParserOk  ("/*  \\*/  */" ,commentParser("/*","*/"), "commentParser")
  checkTokenParserOk  ("\\*\\*"       ,commentParser("\\*","\\*"), "commentParser")
  checkTokenParserOk  ("\\* sdas \\*" ,commentParser("\\*","\\*"), "commentParser")
  checkTokenParserOk  ("\\* \\\\* \\*",commentParser("\\*","\\*"), "commentParser")


  checkTokenParserOk  ("/** **/",commentParser("/*", "*/"), "commentParser")
  checkTokenParserOk  ("/****/" ,commentParser("/*", "*/"), "commentParser")
  checkTokenParserOk  ("/** */" ,commentParser("/*", "*/"), "commentParser")
  checkTokenParserOk  ("/***/"  ,commentParser("/*", "*/"), "commentParser")
  checkTokenParserOk  ("/* \\*/ */"  ,commentParser("/*", "*/"), "commentParser")
  checkTokenParserOk  ("/* *\\/ */"  ,commentParser("/*", "*/"), "commentParser")

  checkTokenParserFail("(*"           ,commentParser("(*","*)"), "commentParser")
  checkTokenParserFail("(* "          ,commentParser("(*","*)"), "commentParser")
  checkTokenParserFail("( * *)"       ,commentParser("(*","*)"), "commentParser")

  checkTokenParserFail("/*  \\*/"     ,commentParser("/*","*/"), "commentParser")
  checkTokenParserFail("\\* \\\\* "   ,commentParser("\\*","\\*"), "commentParser")

  checkTokenParserOk  ("else"       ,keyword("else"),"keyword")

  checkTokenParserOk  (".."        ,dots(),"dots")
  checkTokenParserOk  ("..."       ,dots(),"dots")
  checkTokenParserFail("a..."      ,dots(),"dots")
  checkTokenParserOk  (""  ,empty(),"empty","")
  checkTokenParserOk  ("12",empty(),"empty","")
  checkTokenParserOk  (":"       ,charParser(":"),"char")
  checkTokenParserFail(";"       ,charParser(":"),"char")
  checkTokenParserOk  ("a"       ,charInSetParser(isLetter),"charInSet")
  checkTokenParserOk  ("a"       ,charInSetParser(function (char) TRUE),"charInSet")

  checkTokenParserFail("a"       ,charInSetParser(isDigit ),"charInSet")
  checkTokenParserFail("a"       ,charInSetParser(function (char) FALSE),"charInSet")
  parser12 <- function(p1,p2) { function(stream) { r <- p1(stream); if (r$status=="ok") p2(r$stream) else return(r)  }} 
  checkTokenParserOk  ("   -1233.e+66"         ,parser12(whitespace(),numberScientific()),"numberScientific" ,"-1233.e+66")
  checkTokenParserFail("asd-1233.e+66"         ,parser12(symbolic()      ,numberScientific()),"numberScientific" ,"-1233.e+66")
  checkTokenParserOk  ("asd+1233.e+66"         ,parser12(symbolic()      ,numberScientific()),"numberScientific" ,"+1233.e+66")
  checkTokenParserOk  ("-1233.e+66asd"         ,parser12(numberScientific(),symbolic()),"symbolic" ,"asd")
  checkTokenParserOk("elseasd"               ,parser12(keyword("else"),symbolic()),"symbolic" ,"asd")
  checkTokenParserOk  ("elseasd"               ,symbolic(),"symbolic" ,"elseasd") 
}
test.tokens02 <- function() {
          parseList <- list(
                          whitespace(),
                          string(),whitespace(),
                          numberScientific(),whitespace(),
                          numberNatural(),whitespace(),
                          numberScientific(),whitespace(),
                          numberScientific(),whitespace(),
                          numberInteger(),whitespace(),
                          numberFloat(),whitespace(),
                          numberFloat(),whitespace(),
                          numberScientific(),whitespace(),
                          numberScientific(),whitespace(),
                          numberScientific(),whitespace(),
                          numberScientific(),whitespace(),
                          symbolic(),whitespace(),
                          string(),whitespace(),
                          eofMark()
                  )
   

          for( i in 1:length(parseList) ) {
            
            stream   <- streamParserFromFileName( system.file( "extdata","datInTest01.txt", package = "qmrparser"))
                       
            for( j in 1:i ) {
              cstream  <- parseList[[j]](stream)
              stream   <- cstream$stream                  
            }
            streamParserClose(stream)
            checkEquals("ok",cstream$status,msg=paste(" hasta ",as.character(i)))
          }
}
