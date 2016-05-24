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
  # Check kinds of tokens in PC-AXIS files
    
  name     <- system.file("extdata","datInSFexample6_1.px", package = "qmrparser")
  rule     <- alternation(symbolic(),string(),numberScientific(),
                          charParser("("),charParser(")"),
                          charParser("["),charParser("]"),
                          charParser(","),charParser(";"),charParser("-"),charParser("="),
                          eofMark(),
                          charParser(' '),
                          charParser('\n'),charParser("\r"),charParser("\t"))

  stream   <- streamParserFromFileName(name,encoding="UTF-8")

  cstream  <- rule(stream)
  printCStream(cstream) ; 
  checkEquals("ok",cstream$status,'tokens read')

  while( cstream$status == "ok" &&  cstream$node$value$type != "eofMark" ) { 
    cstream  <- rule(cstream$stream)    
    printCStream(cstream) ;     
    checkEquals("ok",cstream$status,'tokens read')
  }

  if( cstream$status == "fail" ) {     
    print(streamParserNextChar(cstream$stream)) ;  
      print(streamParserPosition(cstream$stream)) ;     
  }
  streamParserClose(cstream$stream)

    
  # step by step PC-AXIS format checking
  checkEquals("ok",pcAxisParser(streamParserFromString('CHARSET="ANSI";'))$status,'CHARSET')

  checkEquals("ok",pcAxisParser(streamParserFromString('AXIS-VERSION="2000"; '))$status,'AXIS-VERSION')

  checkEquals("ok",pcAxisParser(streamParserFromString('DECIMALS=0;'))$status,'DECIMALS')

  checkEquals("ok",pcAxisParser(streamParserFromString('SUBJECT-AREA="Väestö";'))$status,'SUBJECT-AREA')

  checkEquals("ok",pcAxisParser(streamParserFromString('SUBJECT-AREA="Väestö";
  '))$status,'SUBJET-AREA/TITLE/CONTENTS')
      
  checkEquals("ok",pcAxisParser(streamParserFromString('TITLE="Väestö 31.12. muuttujina Sukupuoli, Kunta, Vuosi ja Siviilisääty";'))$status,'TITLE')
   
  checkEquals("ok",pcAxisParser(streamParserFromString('SUBJECT-AREA="Väestö";
  TITLE="Väestö 31.12. muuttujina Sukupuoli, Kunta, Vuosi ja Siviilisääty";
  CONTENTS="Väestö 31.12.";'))$status,'SUBJET-AREA/TITLE/CONTENTS')

  checkEquals("ok",pcAxisParser(streamParserFromString('HEADING="Vuosi","Siviilisääty";'))$status,'HEADING')

  checkEquals("ok",pcAxisParser(streamParserFromString('VALUES("Kunta")="Espoo","Helsinki","Vantaa";'))$status,'VALUES')

  checkEquals("ok",pcAxisParser(streamParserFromString('DATA=57516 43030 100546 57516 43030 100546 91202 67623 158825 91202 67623 158825 ; ' ))$status,'DATA 01')
         
  checkEquals("ok",pcAxisParser(streamParserFromString('DATA=57516 43030 100546 57516 43030 100546 91202 67623 158825 91202 67623 158825  ' ))$status,'DATA 02')
          
  checkEquals("ok",pcAxisParser(streamParserFromString('DATA=
  57516 43030 100546 57516 43030 100546
  144564 85339 229903 144564 85339 229903
  47131 33688 80819 47131 33688 80819
  53821 43375 97196 53821 43375 97196
  151536 86184 237720 151536 86184 237720
  44071 33935 78006 44071 33935 78006
  111337 86405 197742 111337 86405 197742
  296100 171523 467623 296100 171523 467623
  91202 67623 158825 91202 67623 158825;'))$status,'DATA')

  checkEquals("ok",pcAxisParser(streamParserFromString('DATA=
  39669394 19399549 20269845
  13782827 6554619 7228208
     52709   28052   24657
    739409  382664  356745
    800097  406813  393284
   1444239  727941  716298
   3129220 1563613 1565607
   3599227 1789972 1809255
   4525296 2241138 2284158
   4996377 2467192 2529185
   3157049 1549638 1607411
   3442944 1687907 1755037
         0       0       0
  \032'))$status,'DATA')
  
test.filesPcAxis01 <- function() { 

  auxfun00 <- function(stream) {
    cstream  <- pcAxisParser(stream)
    streamParserClose(cstream$stream)
    print(names(cstream))                        
    checkEquals("ok",cstream$status,name)
    
    if( cstream$status != "ok" )  {
      printCStream(cstream)                        
    } else {
      cube <- pcAxisCubeMake(cstream)
      print(names(cube))
    }
  }
        
  auxfun01 <- function(pcaxisFiles,encoding= getOption("encoding"))  {
    for ( name in pcaxisFiles) { 
      
      print("");print("");
      print(paste("File:",name))
      
      name <- system.file("extdata", name, package = "qmrparser")

      stream   <- streamParserFromFileName(name,encoding=encoding)
    
  # print(stream)
      auxfun00(stream)
    }
  }

  auxfun02 <- function(pcaxisFiles,encoding)  {
    for ( name in pcaxisFiles) { 
      
      print("");print("");
      print(paste("File:",name))
      
      name <- system.file("extdata", name, package = "qmrparser")
      
      stream   <- streamParserFromString(iconv(readLines(name,encoding),"UTF-8"))
  # print(stream)
      auxfun00(stream)      
    }
  }
  
  ## from sample files:
  pcaxisFilesExamples <- list(
                      "datInSFexample6_1.px", "datInSFexample6_2.px" ,
                      "datInSFexample6_3.px", "datInSFexample6_4.px" , 
                      "datInSFexampleA_5.px", "datInSFexample6_5.px" )
  auxfun01(pcaxisFilesExamples,"UTF-8")
  


} # end function
