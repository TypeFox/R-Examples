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
test.rules01 <- function() {
  # alternation
  checkSentParserOk  (
  "'aaa 123'",
  alternation(string(action=actionString),numberNatural(action=action)),
  "alternation")

  checkSentParserOk  (
  "1234",
  alternation(string(action=actionString),numberNatural(action=action)),
  "alternation")
  # option
  checkSentParserOk  (
  "'aaa1234'",
  option(string(action=actionString)),
  "option")

  checkSentParserOk  (
  "1234",
  option(string(action=actionString)),
  "option","empty")

  checkSentParserOk  (
  "1234",
  option(string(action=actionString),action=function(s) if( is.list(s) ) list(type="option",value=s$value) else list(type="option",value=s)),
  "option","")

  checkSentParserOk  (
  "'1234'",
  option(string(action=actionString),action=function(s) if( is.list(s) ) list(type="option",value=s$value) else list(type="option",value=s)),
  "option")
  # concatenation
  checkSentParserOk  (
  " 123",
  concatenation(whitespace(action=action),numberNatural(action=action)),
  "concatenation")

  checkSentParserOk  (
  " 'abs  21'",
  concatenation(whitespace(action=action),string(action=actionString)),
  "concatenation")

  checkSentParserOk  (
  " a123",
  concatenation(whitespace(action=action),symbolic(action=action)),
  "concatenation")

  checkSentParserOk  (
  " 123 abcd",
  concatenation(whitespace(action=action),numberNatural(action=action),whitespace(action=action),symbolic(action=action)),
  "concatenation")

  checkSentParserOk  (
  " 123abcd",
  concatenation(whitespace(action=action),numberNatural(action=action),whitespace(action=action),symbolic(action=action)),
  "concatenation")

  checkSentParserOk  (
  " 'aa' 123",
  concatenation(whitespace(action=action),string(action=actionString),whitespace(action=action),numberNatural(action=action)),
  "concatenation")

  checkSentParserFail(
  "1 'aa' 123",
  concatenation(whitespace(action=action),string(action=actionString),whitespace(action=action),numberNatural(action=action)),
  "concatenation")
  # repetition1N
  checkSentParserOk  (
  "'1234'",
  repetition1N(string(action=actionString)),
  "repetition1N")

  checkSentParserOk  (
  "'1234''1234'",
  repetition1N(string(action=actionString)),
  "repetition1N")

  checkSentParserOk  (
  "'1234''1234''1234'",
  repetition1N(string(action=actionString)),
  "repetition1N")
  # repetition0N

  checkSentParserOk  (
  "1234",
  repetition0N(string(action=actionString)),
  "repetition0N","empty")

  checkSentParserOk  (
  "'1234'",
  repetition0N(string(action=actionString)),
  "repetition0N",
  "repetition1Nlist(\"'1234'\")"
  )
  checkSentParserOk  (
  "'1234''1234'",
  repetition0N(string(action=actionString)),
  "repetition0N",
  "repetition1Nlist(\"'1234'\", \"'1234'\")"
  )
  checkSentParserOk  (
  "'1234''1234''1234'",
  repetition0N(string(action=actionString)),
  "repetition0N",
  "repetition1Nlist(\"'1234'\", \"'1234'\", \"'1234'\")"
  )

  checkSentParserOk  (
  "'1234''1234''1234'",
  repetition0N(string(action=actionString),action= function(s) if( is.list(s) ) list(type="repetition0N",value=s$value) else list(type="repetition0N",value=s)),
  "repetition0N",
  )

  checkSentParserOk  (
  "1234'1234''1234'",
  repetition0N(string(action=actionString),action= function(s) if( is.list(s) ) list(type="repetition0N",value=s$value) else list(type="repetition0N",value=s)),
  "repetition0N",""
  )
}
test.rules02 <- function() {
  # combinations

  checkSentParserOk  (
  " 'aa'",
  concatenation(whitespace(action=action),alternation(string(action=actionString),numberNatural(action=action),action=action)),
  "concatenation")

  checkSentParserOk  (
  " 123",
  concatenation(whitespace(action=action),alternation(string(action=actionString),numberNatural(action=action),action=action)),
  "concatenation")

  checkSentParserOk  (
  " 123  ",
  concatenation(whitespace(action=action),alternation(string(action=actionString),numberNatural(action=action),action=action),whitespace(action=action)),
  "concatenation")
          parser <- concatenation(
                          whitespace(),
                          string(),whitespace(),
                          repetition1N(concatenation(numberScientific(),whitespace())),
                          symbolic(),whitespace(),
                          string(),whitespace(),
                          eofMark()
                  )

          stream   <- streamParserFromFileName(system.file("extdata","datInTest01.txt", package = "qmrparser"))

          cstream  <-  parser(stream)
          streamParserClose(cstream$stream)
          checkEquals("ok",cstream$status)

}
