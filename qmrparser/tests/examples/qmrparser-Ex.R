pkgname <- "qmrparser"
source(file.path(R.home("share"), "R", "examples-header.R"))
options(warn = 1)
library('qmrparser')

base::assign(".oldSearch", base::search(), pos = 'CheckExEnv')
cleanEx()
nameEx("alternation")
### * alternation

flush(stderr()); flush(stdout())

### Name: alternation
### Title: Alternative phrases
### Aliases: alternation
### Keywords: parser combinator

### ** Examples



# ok
stream  <- streamParserFromString("123 Hello world")
( alternation(numberNatural(),symbolic())(stream) )[c("status","node")]


# fail
stream  <- streamParserFromString("123 Hello world")
( alternation(string(),symbolic())(stream) )[c("status","node")]





cleanEx()
nameEx("charInSetParser")
### * charInSetParser

flush(stderr()); flush(stdout())

### Name: charInSetParser
### Title: Single character, belonging to a given set, token
### Aliases: charInSetParser
### Keywords: token

### ** Examples


# fail
stream  <- streamParserFromString("H")
( charInSetParser(isDigit)(stream) )[c("status","node")]

# ok 
stream  <- streamParserFromString("a")
( charInSetParser(isLetter)(stream) )[c("status","node")]





cleanEx()
nameEx("charParser")
### * charParser

flush(stderr()); flush(stdout())

### Name: charParser
### Title: Specific single character token.
### Aliases: charParser
### Keywords: token

### ** Examples


# fail
stream  <- streamParserFromString("H")
( charParser("a")(stream) )[c("status","node")]

# ok 
stream  <- streamParserFromString("a")
( charParser("a")(stream) )[c("status","node")]

# ok 
( charParser("\U00B6")(streamParserFromString("\U00B6")) )[c("status","node")]




cleanEx()
nameEx("commentParser")
### * commentParser

flush(stderr()); flush(stdout())

### Name: commentParser
### Title: Comment token.
### Aliases: commentParser
### Keywords: token

### ** Examples


# fail
stream  <- streamParserFromString("123")
( commentParser("(*","*)")(stream) )[c("status","node")]

# ok 
stream  <- streamParserFromString("(*123*)")
( commentParser("(*","*)")(stream) )[c("status","node")]




cleanEx()
nameEx("concatenation")
### * concatenation

flush(stderr()); flush(stdout())

### Name: concatenation
### Title: One phrase then another
### Aliases: concatenation
### Keywords: parser combinator

### ** Examples


# ok
stream  <- streamParserFromString("123Hello world")
( concatenation(numberNatural(),symbolic())(stream) )[c("status","node")]


# fail
stream  <- streamParserFromString("123 Hello world")
( concatenation(string(),symbolic())(stream) )[c("status","node")]




cleanEx()
nameEx("dots")
### * dots

flush(stderr()); flush(stdout())

### Name: dots
### Title: Dots sequence token.
### Aliases: dots
### Keywords: token

### ** Examples


# fail
stream  <- streamParserFromString("Hello world")
( dots()(stream) )[c("status","node")]

# ok 
stream  <- streamParserFromString("..")
( dots()(stream) )[c("status","node")]




cleanEx()
nameEx("empty")
### * empty

flush(stderr()); flush(stdout())

### Name: empty
### Title: Empty token
### Aliases: empty
### Keywords: token

### ** Examples


# ok
stream  <- streamParserFromString("Hello world")
( empty()(stream) )[c("status","node")]

# ok 
stream  <- streamParserFromString("")
( empty()(stream) )[c("status","node")]




cleanEx()
nameEx("eofMark")
### * eofMark

flush(stderr()); flush(stdout())

### Name: eofMark
### Title: End of file token
### Aliases: eofMark
### Keywords: token

### ** Examples


# fail
stream  <- streamParserFromString("Hello world")
( eofMark()(stream) )[c("status","node")]

# ok 
stream  <- streamParserFromString("")
( eofMark()(stream) )[c("status","node")]




cleanEx()
nameEx("isDigit")
### * isDigit

flush(stderr()); flush(stdout())

### Name: isDigit
### Title: Is it a digit?
### Aliases: isDigit
### Keywords: set of character

### ** Examples

isDigit('9')
isDigit('a')



cleanEx()
nameEx("isHex")
### * isHex

flush(stderr()); flush(stdout())

### Name: isHex
### Title: Is it an hexadecimal digit?
### Aliases: isHex
### Keywords: set of character

### ** Examples

isHex('+')
isHex('A')
isHex('a')
isHex('9')



cleanEx()
nameEx("isLetter")
### * isLetter

flush(stderr()); flush(stdout())

### Name: isLetter
### Title: Is it a letter?
### Aliases: isLetter
### Keywords: set of character

### ** Examples

isLetter('A')
isLetter('a')
isLetter('9')



cleanEx()
nameEx("isLowercase")
### * isLowercase

flush(stderr()); flush(stdout())

### Name: isLowercase
### Title: Is it a lower case?
### Aliases: isLowercase
### Keywords: set of character

### ** Examples

isLowercase('A')
isLowercase('a')
isLowercase('9')



cleanEx()
nameEx("isNewline")
### * isNewline

flush(stderr()); flush(stdout())

### Name: isNewline
### Title: Is it a new line character?
### Aliases: isNewline
### Keywords: set of character

### ** Examples

isNewline(' ')
isNewline('\n')



cleanEx()
nameEx("isSymbol")
### * isSymbol

flush(stderr()); flush(stdout())

### Name: isSymbol
### Title: Is it a symbol?
### Aliases: isSymbol
### Keywords: set of character

### ** Examples

isSymbol('+')
isSymbol('A')
isSymbol('a')
isSymbol('9')



cleanEx()
nameEx("isUppercase")
### * isUppercase

flush(stderr()); flush(stdout())

### Name: isUppercase
### Title: Is it an upper case?
### Aliases: isUppercase
### Keywords: set of character

### ** Examples

isUppercase('A')
isUppercase('a')
isUppercase('9')



cleanEx()
nameEx("isWhitespace")
### * isWhitespace

flush(stderr()); flush(stdout())

### Name: isWhitespace
### Title: Is it a white space?
### Aliases: isWhitespace
### Keywords: set of character

### ** Examples

isWhitespace(' ')
isWhitespace('\n')
isWhitespace('a')



cleanEx()
nameEx("keyword")
### * keyword

flush(stderr()); flush(stdout())

### Name: keyword
### Title: Arbitrary given token.
### Aliases: keyword
### Keywords: token

### ** Examples


# fail
stream  <- streamParserFromString("Hello world")
( keyword("world")(stream) )[c("status","node")]

# ok 
stream  <- streamParserFromString("world")
( keyword("world")(stream) )[c("status","node")]




cleanEx()
nameEx("numberFloat")
### * numberFloat

flush(stderr()); flush(stdout())

### Name: numberFloat
### Title: Floating-point number token.
### Aliases: numberFloat
### Keywords: token

### ** Examples


# fail
stream  <- streamParserFromString("Hello world")
( numberFloat()(stream) )[c("status","node")]

# ok 
stream  <- streamParserFromString("-456.74")
( numberFloat()(stream) )[c("status","node")]




cleanEx()
nameEx("numberInteger")
### * numberInteger

flush(stderr()); flush(stdout())

### Name: numberInteger
### Title: Integer number token.
### Aliases: numberInteger
### Keywords: token

### ** Examples


# fail
stream  <- streamParserFromString("Hello world")
( numberInteger()(stream) )[c("status","node")]

# ok 
stream  <- streamParserFromString("-1234")
( numberInteger()(stream) )[c("status","node")]




cleanEx()
nameEx("numberNatural")
### * numberNatural

flush(stderr()); flush(stdout())

### Name: numberNatural
### Title: Natural number token.
### Aliases: numberNatural
### Keywords: token

### ** Examples


# fail
stream  <- streamParserFromString("Hello world")
( numberNatural()(stream) )[c("status","node")]

# ok 
stream  <- streamParserFromString("123")
( numberNatural()(stream) )[c("status","node")]




cleanEx()
nameEx("numberScientific")
### * numberScientific

flush(stderr()); flush(stdout())

### Name: numberScientific
### Title: Number in scientific notation token.
### Aliases: numberScientific
### Keywords: token

### ** Examples


# fail
stream  <- streamParserFromString("Hello world")
( numberScientific()(stream) )[c("status","node")]

# ok 
stream  <- streamParserFromString("-1234e12")
( numberScientific()(stream) )[c("status","node")]




cleanEx()
nameEx("option")
### * option

flush(stderr()); flush(stdout())

### Name: option
### Title: Optional parser
### Aliases: option
### Keywords: parser combinator

### ** Examples


# ok
stream  <- streamParserFromString("123 Hello world")
( option(numberNatural())(stream) )[c("status","node")]


# ok
stream  <- streamParserFromString("123 Hello world")
( option(string())(stream) )[c("status","node")]




cleanEx()
nameEx("pcAxisCubeMake")
### * pcAxisCubeMake

flush(stderr()); flush(stdout())

### Name: pcAxisCubeMake
### Title: Creates PC-AXIS cube
### Aliases: pcAxisCubeMake
### Keywords: PC-AXIS

### ** Examples


  ## Not run: 
##D     ## significant time reductions may be achieve by doing:
##D     library("compiler")
##D     enableJIT(level=3)
##D   
## End(Not run)
  
  name     <- system.file("extdata","datInSFexample6_1.px", package = "qmrparser")
  
  stream   <- streamParserFromFileName(name,encoding="UTF-8")
  
  cstream  <-  pcAxisParser(stream)
  if ( cstream$status == 'ok' ) {
    cube <- pcAxisCubeMake(cstream)
    
    ## Variables
    print(cube$pxCubeVariable)
    
    ## Data
    print(cube$pxCubeData)

  }
  
  ## Not run: 
##D       #
##D       # Error messages like
##D       #                " ... invalid multibyte string ... "
##D       # or warnings
##D       #                " input string ...  is invalid in this locale"
##D       #
##D       # For example, in Linux the error generated by this code:
##D        name     <-     "http://www.ine.es/pcaxisdl//t20/e245/p04/a2009/l0/00000008.px" 
##D       stream   <- streamParserFromString( readLines( name ) )    
##D       cstream  <- pcAxisParser(stream)
##D       if ( cstream$status == 'ok' )  cube <- pcAxisCubeMake(cstream)
##D       #
##D       # is caused by files with a non-readable 'encoding'.
##D       # In the case where it could be read, there may also be problems 
##D       # with string-handling functions, due to multibyte characters. 
##D       # In Windows, according to \code{link{Sys.getlocale}()},
##D       # file may be read but accents, ñ, ... may not be correctly recognised.
##D       #
##D       #
##D       # There are, at least, the following options:
##D       #  - File conversion to utf-8, from the OS, with
##D       # "iconv - Convert encoding of given files from one encoding to another"
##D       #
##D       #  - File conversion in R:
##D       name    <- "http://www.ine.es/pcaxisdl//t20/e245/p04/a2009/l0/00000008.px" 
##D       stream   <- streamParserFromString( iconv( readLines( name ), "IBM850", "UTF-8") )
##D       cstream  <- pcAxisParser(stream)
##D       if ( cstream$status == 'ok' )  cube <- pcAxisCubeMake(cstream)
##D       #
##D       # In the latter case, latin1 would also work, but accents, ñ, ... would not be correctly read.
##D       #
##D       #  - Making the assumption that the file does not contain multibyte characters:
##D       #
##D       localeOld <- Sys.getlocale("LC_CTYPE")
##D       Sys.setlocale(category = "LC_CTYPE", locale = "C")
##D       #
##D       name     <-
##D         "http://www.ine.es/pcaxisdl//t20/e245/p04/a2009/l0/00000008.px" 
##D       stream   <- streamParserFromString( readLines( name ) )
##D       cstream  <- pcAxisParser(stream)
##D       if ( cstream$status == 'ok' )  cube <- pcAxisCubeMake(cstream)
##D       #
##D       Sys.setlocale(category = "LC_CTYPE", locale = localeOld)
##D       #
##D       # However, some characters will not be correctly read (accents, ñ, ...)
##D 
##D     
## End(Not run)
    



cleanEx()
nameEx("pcAxisCubeToCSV")
### * pcAxisCubeToCSV

flush(stderr()); flush(stdout())

### Name: pcAxisCubeToCSV
### Title: Exports a PC-AXIS cube into CSV in several files.
### Aliases: pcAxisCubeToCSV
### Keywords: PC-AXIS

### ** Examples


  name     <- system.file("extdata","datInSFexample6_1.px", package = "qmrparser")
  stream   <- streamParserFromFileName(name,encoding="UTF-8")
  cstream  <-  pcAxisParser(stream)
  if ( cstream$status == 'ok' ) {
    cube <- pcAxisCubeMake(cstream)
    
    pcAxisCubeToCSV(prefix="datInSFexample6_1",pcAxisCube=cube)     

  }



cleanEx()
nameEx("pcAxisParser")
### * pcAxisParser

flush(stderr()); flush(stdout())

### Name: pcAxisParser
### Title: Parser for PC-AXIS format files
### Aliases: pcAxisParser
### Keywords: PC-AXIS

### ** Examples


  ## Not run: 
##D     ## significant time reductions may be achieve by doing:
##D     library("compiler")
##D     enableJIT(level=3)
##D   
## End(Not run)

  name     <- system.file("extdata","datInSFexample6_1.px", package = "qmrparser")
  stream   <- streamParserFromFileName(name,encoding="UTF-8")
  cstream  <-  pcAxisParser(stream)
  if ( cstream$status == 'ok' ) {
    
    ## HEADING 
    print(Filter(function(e) e$keyword=="HEADING",cstream$node)[[1]] $ruleRight$value)  
  
    ## STUB
    print(Filter(function(e) e$keyword=="STUB",cstream$node)[[1]] $ruleRight$value)  
  
    ## DATA
    print(Filter(function(e) e$keyword=="DATA",cstream$node)[[1]] $ruleRight$value)
    
  }

  ## Not run: 
##D       #
##D       # Error messages like
##D       #                " ... invalid multibyte string ... "
##D       # or warnings
##D       #                " input string ...  is invalid in this locale"
##D       #
##D       # For example, in Linux the error generated by this code:
##D        name     <-     "http://www.ine.es/pcaxisdl//t20/e245/p04/a2009/l0/00000008.px" 
##D       stream   <- streamParserFromString( readLines( name ) )    
##D       cstream  <- pcAxisParser(stream)
##D       if ( cstream$status == 'ok' )  cube <- pcAxisCubeMake(cstream)
##D       #
##D       # is caused by files with a non-readable 'encoding'.
##D       # In the case where it could be read, there may also be problems 
##D       # with string-handling functions, due to multibyte characters. 
##D       # In Windows, according to \code{link{Sys.getlocale}()},
##D       # file may be read but accents, ñ, ... may not be correctly recognised.
##D       #
##D       #
##D       # There are, at least, the following options:
##D       #  - File conversion to utf-8, from the OS, with
##D       # "iconv - Convert encoding of given files from one encoding to another"
##D       #
##D       #  - File conversion in R:
##D       name    <- "http://www.ine.es/pcaxisdl//t20/e245/p04/a2009/l0/00000008.px" 
##D       stream   <- streamParserFromString( iconv( readLines( name ), "IBM850", "UTF-8") )
##D       cstream  <- pcAxisParser(stream)
##D       if ( cstream$status == 'ok' )  cube <- pcAxisCubeMake(cstream)
##D       #
##D       # In the latter case, latin1 would also work, but accents, ñ, ... would not be correctly read.
##D       #
##D       #  - Making the assumption that the file does not contain multibyte characters:
##D       #
##D       localeOld <- Sys.getlocale("LC_CTYPE")
##D       Sys.setlocale(category = "LC_CTYPE", locale = "C")
##D       #
##D       name     <-
##D         "http://www.ine.es/pcaxisdl//t20/e245/p04/a2009/l0/00000008.px" 
##D       stream   <- streamParserFromString( readLines( name ) )
##D       cstream  <- pcAxisParser(stream)
##D       if ( cstream$status == 'ok' )  cube <- pcAxisCubeMake(cstream)
##D       #
##D       Sys.setlocale(category = "LC_CTYPE", locale = localeOld)
##D       #
##D       # However, some characters will not be correctly read (accents, ñ, ...)
##D 
##D     
## End(Not run)




cleanEx()
nameEx("repetition0N")
### * repetition0N

flush(stderr()); flush(stdout())

### Name: repetition0N
### Title: Repeats one parser
### Aliases: repetition0N
### Keywords: parser combinator

### ** Examples


# ok
stream  <- streamParserFromString("Hello world")
( repetition0N(symbolic())(stream) )[c("status","node")]


# ok
stream  <- streamParserFromString("123 Hello world")
( repetition0N(symbolic())(stream) )[c("status","node")]




cleanEx()
nameEx("repetition1N")
### * repetition1N

flush(stderr()); flush(stdout())

### Name: repetition1N
### Title: Repeats a parser, at least once.
### Aliases: repetition1N
### Keywords: parser combinator

### ** Examples


# ok
stream  <- streamParserFromString("Hello world")
( repetition1N(symbolic())(stream) )[c("status","node")]


# fail
stream  <- streamParserFromString("123 Hello world")
( repetition1N(symbolic())(stream) )[c("status","node")]




cleanEx()
nameEx("separator")
### * separator

flush(stderr()); flush(stdout())

### Name: separator
### Title: Generic word separator token.
### Aliases: separator
### Keywords: token

### ** Examples


# ok
stream  <- streamParserFromString("; Hello world")
( separator()(stream) )[c("status","node")]

# ok 
stream  <- streamParserFromString(" ")
( separator()(stream) )[c("status","node")]


# fail
stream  <- streamParserFromString("Hello world")
( separator()(stream) )[c("status","node")]

# fail 
stream  <- streamParserFromString("")
( separator()(stream) )[c("status","node")]




cleanEx()
nameEx("streamParser")
### * streamParser

flush(stderr()); flush(stdout())

### Name: streamParser
### Title: Generic interface for character processing, allowing forward and
###   backwards translation.
### Aliases: streamParserNextChar streamParserNextChar
###   streamParserNextCharSeq streamParserPosition streamParserClose
### Keywords: streamParser

### ** Examples


stream<- streamParserFromString("Hello world")

cstream <- streamParserNextChar(stream)

while( cstream$status == "ok" ) {
    print(streamParserPosition(cstream$stream))
    print(cstream$char)
    cstream <- streamParserNextCharSeq(cstream$stream)
}

streamParserClose(stream)




cleanEx()
nameEx("streamParserFromFileName")
### * streamParserFromFileName

flush(stderr()); flush(stdout())

### Name: streamParserFromFileName
### Title: Creates a streamParser from a file name
### Aliases: streamParserFromFileName
### Keywords: streamParser

### ** Examples

  name    <- system.file("extdata","datInTest01.txt", package = "qmrparser")
  
  stream  <- streamParserFromFileName(name)
  
  cstream <- streamParserNextChar(stream)
  
  while( cstream$status == "ok" ) {
    print(streamParserPosition(cstream$stream))
    print(cstream$char)
    cstream <- streamParserNextCharSeq(cstream$stream)
  }
  
  streamParserClose(stream)
  



cleanEx()
nameEx("streamParserFromString")
### * streamParserFromString

flush(stderr()); flush(stdout())

### Name: streamParserFromString
### Title: Creates a streamParser from a string
### Aliases: streamParserFromString
### Keywords: streamParser

### ** Examples

# reads one character
streamParserNextChar(streamParserFromString("\U00B6"))

# reads a string
stream  <- streamParserFromString("Hello world")

cstream <- streamParserNextChar(stream)

while( cstream$status == "ok" ) {
    print(streamParserPosition(cstream$stream))
    print(cstream$char)
    cstream <- streamParserNextCharSeq(cstream$stream)

streamParserClose(stream)
}




cleanEx()
nameEx("string")
### * string

flush(stderr()); flush(stdout())

### Name: string
### Title: Token string
### Aliases: string
### Keywords: token

### ** Examples


# fail
stream  <- streamParserFromString("Hello world")
( string()(stream) )[c("status","node")]

# ok 
stream  <- streamParserFromString("'Hello world'")
( string()(stream) )[c("status","node")]




cleanEx()
nameEx("symbolic")
### * symbolic

flush(stderr()); flush(stdout())

### Name: symbolic
### Title: Alphanumeric token.
### Aliases: symbolic
### Keywords: token

### ** Examples


# fail
stream  <- streamParserFromString("123")
( symbolic()(stream) )[c("status","node")]

# ok 
stream  <- streamParserFromString("abc123_2")
( symbolic()(stream) )[c("status","node")]




cleanEx()
nameEx("whitespace")
### * whitespace

flush(stderr()); flush(stdout())

### Name: whitespace
### Title: White sequence token.
### Aliases: whitespace
### Keywords: token

### ** Examples


# ok
stream  <- streamParserFromString("Hello world")
( whitespace()(stream) )[c("status","node")]

# ok
stream  <- streamParserFromString(" Hello world")
( whitespace()(stream) )[c("status","node")]

# ok 
stream  <- streamParserFromString("")
(  whitespace()(stream) )[c("status","node")]




### * <FOOTER>
###
options(digits = 7L)
base::cat("Time elapsed: ", proc.time() - base::get("ptime", pos = 'CheckExEnv'),"\n")
grDevices::dev.off()
###
### Local variables: ***
### mode: outline-minor ***
### outline-regexp: "\\(> \\)?### [*]+" ***
### End: ***
quit('no')
