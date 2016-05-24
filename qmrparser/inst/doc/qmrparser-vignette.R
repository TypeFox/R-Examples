### R code from vignette source 'qmrparser-vignette.Rnw'
### Encoding: UTF-8

###################################################
### code chunk number 1: vignette
###################################################
list_expression  <- function() 
concatenation(
concatenation(additive_expression(),ws(),keyword(';'),
action=function(s) print(exprToString(s[[1]]))),
repetition0N(
concatenation(additive_expression(),ws(),keyword(';'),
action=function(s) print(exprToString(s[[1]])))),
action=function(s) NULL)


###################################################
### code chunk number 2: vignette
###################################################
library(qmrparser)

additive_expression  <- function()  concatenation(ws(),multiplicative_expression(),
option( concatenation(ws(),
alternation(
keyword('+',action=function(s) s),
keyword('-',action=function(s) s),
action=function(s) s),
ws(),additive_expression(),
action=function(s) list(type='noempty',value=s[c(2,4)]))),
action=function(s) {if(s[[3]]$value$type=='empty') s[[2]] else
list(fun=s[[3]]$value$value[[1]],par1=s[[2]],par2=s[[3]]$value$value[[2]])
});


###################################################
### code chunk number 3: vignette
###################################################
multiplicative_expression  <- function()  
concatenation(power_expression(),
option(
concatenation(ws(),
alternation(
keyword('*',action=function(s) s),
keyword('/',action=function(s) s),
action=function(s) s),
ws(),multiplicative_expression(),
action=function(s) {list(type='noempty',value=s[c(2,4)])})),
action=function(s) {
if(s[[2]]$value$type=='empty') s[[1]] else
list(fun=s[[2]]$value$value[[1]],par1=s[[1]],par2=s[[2]]$value$value[[2]])});


###################################################
### code chunk number 4: vignette
###################################################
power_expression  <- function()  
concatenation(primary(),
option(
concatenation(ws(),keyword('**'),ws(),power_expression(),
action=function(s) list(type='noempty',value=s[[4]]))),
action=function(s){if(s[[2]]$value$type=='empty') s[[1]] else list(fun="^",par1=s[[1]],par2=s[[2]]$value$value)});


###################################################
### code chunk number 5: vignette
###################################################
primary  <- function()  alternation(
concatenation(charParser('('),ws(),additive_expression(),ws(),charParser(')'),action=function(s)  s[[3]]),

concatenation(charParser('-'),ws(),primary(),
action=function(s) list(fun="U-",par1=s[[3]])),

concatenation(FUN(action=function(s) s), ws(), charParser('('), ws(), additive_expression(), ws(), charParser(')'),
action=function(s) list(fun=s[[1]],par1=s[[5]])),

NUMBER  (action=function(s) list(fun="NUM", par1=s)),

VARIABLE(action=function(s) list(fun="VAR", par1=s)),
action=function(s) s);


###################################################
### code chunk number 6: vignette
###################################################
NUMBER    <- function(...)  numberScientific(...);

VARIABLE  <- function(...)  symbolic(...);

FUN       <- function(...)  symbolic(...);

ws        <- function()     whitespace();


###################################################
### code chunk number 7: vignette
###################################################
exprToString <- function(expr) 
if ( !is.list(expr) ) as.character(expr) else
paste("(",paste(sapply(expr,exprToString,USE.NAMES=FALSE),collapse=" "),")")


###################################################
### code chunk number 8: vignette
###################################################
print("Infix to Prefix Examples")
 
invisible( list_expression()(streamParserFromString(" 8  ;")) )

invisible( list_expression()(streamParserFromString("8 +4;")) )

invisible( list_expression()(streamParserFromString("8/2 ;")) )

invisible( list_expression()(streamParserFromString("8*2 ;")) )

invisible( list_expression()(streamParserFromString("2*3 + 4*5;")) )

invisible( list_expression()(streamParserFromString("sqrt( 16) ;")) )

invisible( list_expression()(streamParserFromString("sin(3.1415) ;")) )

invisible( list_expression()(streamParserFromString("sin(3.14* (2*2+3+1)/2 ) ** 8;")) )

invisible( list_expression()(streamParserFromString("sqrt(16)**2+sin(3)-sin(3);")) )

invisible( list_expression()(streamParserFromString("sqrt(16)**2+sin(3)-sin(3)*2;") ) )
           


###################################################
### code chunk number 9: vignette
###################################################
exprToNumber <- function(expr)  
switch(expr[[1]],
'NUM'= as.numeric(expr[[2]]),
'VAR' =as.numeric(get(expr[[2]])),
'U-'=-exprToNumber(expr[[2]]),
do.call(expr[[1]],unname(lapply(expr[-1],exprToNumber)))
)	 


###################################################
### code chunk number 10: vignette
###################################################
list_expression  <- function() 
concatenation(
concatenation(additive_expression(),ws(),keyword(';'),
action=function(s) print(exprToString(exprToNumber(s[[1]])))),
repetition0N(
concatenation(additive_expression(),ws(),keyword(';'),
action=function(s) print(exprToString(exprToNumber(s[[1]]))))),
action=function(s) NULL);


###################################################
### code chunk number 11: vignette
###################################################
print("Calculator")

invisible( list_expression()(streamParserFromString(" 8  ;")) )
          
invisible( list_expression()(streamParserFromString("8 +4;")) )

invisible( list_expression()(streamParserFromString("8/2 ;")) )
          
invisible( list_expression()(streamParserFromString("8*2 ;")) )
          
invisible( list_expression()(streamParserFromString("2*3 + 4*5;")) )

invisible( list_expression()(streamParserFromString("sqrt( 16) ;")) )

invisible( list_expression()(streamParserFromString("sin(3.1415) ;")) )

invisible( list_expression()(streamParserFromString("sin(3.14* (2*2+3+1)/2 ) ** 8;")) )

invisible( list_expression()(streamParserFromString("sqrt(16)**2+sin(3)-sin(3);")) )

invisible( list_expression()(streamParserFromString("sqrt(16)**2+sin(3)-sin(3)*2;")) )



###################################################
### code chunk number 12: vignette
###################################################
exprDeriv <- function(expr,var)  
switch(expr[[1]],
'NUM'= list("NUM", "0"),
'VAR'  = if( expr[[2]] == var ) list("NUM" ,"1") else list("NUM", "0"),
"+"=,"-"= list(expr[[1]],exprDeriv(expr[[2]],var),exprDeriv(expr[[3]],var)),
"*"    =list("+",
              list("*",expr[[2]],exprDeriv(expr[[3]],var)),
              list("*",expr[[3]],exprDeriv(expr[[2]],var))
	 ),
"/"    =list("*",
             list("-",
              list("*",expr[[3]],exprDeriv(expr[[2]],var)),
              list("*",expr[[2]],exprDeriv(expr[[3]],var))
	      
	     ),	 
	     list("**",expr[[3]],"2")
	 ),
"sin"=list("*",exprDeriv(expr[[2]],var),list("cos", expr[[2]])),
	 list(paste("Diff",var,sep="_"),expr)
)	 


###################################################
### code chunk number 13: vignette
###################################################
list_expression  <- function() 
concatenation(
concatenation(additive_expression(),ws(),keyword(';'),
action=function(s) print(exprToString(exprDeriv(s[[1]],"x")))),
repetition0N(
concatenation(additive_expression(),ws(),keyword(';'),
action=function(s) print(exprToString(exprDeriv(s[[1]],"x"))))),
action=function(s) NULL);


###################################################
### code chunk number 14: vignette
###################################################
print("Differentiation")

invisible( list_expression()(streamParserFromString(" 8  ;")) )

invisible( list_expression()(streamParserFromString(" x  ;")) )
          
invisible( list_expression()(streamParserFromString("8 +x;")) )
          
invisible( list_expression()(streamParserFromString("x/2 ;")) )

invisible( list_expression()(streamParserFromString("8*x ;")) )

invisible( list_expression()(streamParserFromString("2*x + 4*x;")) )

invisible( list_expression()(streamParserFromString("1+sqrt( x) ;")) )

invisible( list_expression()(streamParserFromString("sin(x) ;")) )

invisible( list_expression()(streamParserFromString("sin(x* (2*2+x+1)/2 ) ** 8;")) )



###################################################
### code chunk number 15: vignette
###################################################
gramatica <- function() 
concatenation( 
repetition1N(
concatenation(ebnfRule(),whitespace(),action=function(s) s[[1]]),
action=function(s) s),
eofMark(error=function(p) errorFun(p,h=NULL,type="eofMark")),
action=function(s) unlist(s[[1]]) ) 


###################################################
### code chunk number 16: vignette
###################################################
ebnfRule <- function()
concatenation(

whitespace(),

symbolic(charFirst=isLetter,charRest=function(ch) isLetter(ch) || isDigit(ch) || ch == "_",action=function(s) s),

whitespace(),charParser("="),

whitespace(),ebnfDefinition(),whitespace(),charParser(';'),whitespace(),

action=function(s) paste(s[[2]]," <- function() ", s[[6]]))


###################################################
### code chunk number 17: ebnfDefinition
###################################################
ebnfDefinition <- function() alternation(
# several alternatives
ebnfAlternation(),
# No alternatives
ebnfNonAlternation(),
action=function(s) s)


###################################################
### code chunk number 18: ebnfAlternation
###################################################
ebnfAlternation <- function() 
concatenation(
ebnfNonAlternation(),

repetition1N(
concatenation(whitespace(),charParser("|"),whitespace(),ebnfNonAlternation(),action=function(s) s[[4]]),action=function(s) s),
action=function(s) paste("alternation(",paste(s[[1]],",",paste(unlist(s[[2]]), collapse=","), sep=""),")",sep=""))


###################################################
### code chunk number 19: vignette
###################################################
ebnfConcatenation <- function()
option(
concatenation(
whitespace(),charParser(","),whitespace(),
ebnfNonAlternation(),
action=function(s) list(type="noempty",value=s[[4]])))


###################################################
### code chunk number 20: ebnfNonAlternation-string (eval = FALSE)
###################################################
## # string
## concatenation(
## string(action=function(s) paste("keyword('",s,"')",sep="")),
## ebnfConcatenation(),
## action=function (s) 
## if(s[[2]]$value$type=="empty") s[[1]]
## else paste("concatenation(",s[[1]],",",s[[2]]$value$value,")",sep=""))


###################################################
### code chunk number 21: ebnfNonAlternation-special (eval = FALSE)
###################################################
## # special sequence
## concatenation(
## ebnfSpecialSequence(),
## ebnfConcatenation(),
## action=function (s) 
## if(s[[2]]$value$type=="empty") s[[1]]
## else paste("concatenation(",s[[1]],",",s[[2]]$value$value,")",sep=""))


###################################################
### code chunk number 22: ebnfNonAlternation-rule (eval = FALSE)
###################################################
## # rule call
## concatenation(
## symbolic(charFirst=isLetter,charRest=function(ch) isLetter(ch) || isDigit(ch) || ch == "_",action=function(s) paste(s,"()",sep="")),
## ebnfConcatenation(),
## action=function (s) 
## if(s[[2]]$value$type=="empty") s[[1]]
## else paste("concatenation(",s[[1]],",",s[[2]]$value$value,")",sep=""))


###################################################
### code chunk number 23: ebnfNonAlternation-grouping (eval = FALSE)
###################################################
## # grouping
## concatenation(
## whitespace(),charParser("("),whitespace(),
## ebnfDefinition(),
## whitespace(),charParser(")"),
## ebnfConcatenation(),
## action=function (s) 
## if(s[[7]]$value$type=="empty") s[[4]]
## else paste("concatenation(",s[[4]],",",s[[7]]$value$value,")",sep=""))


###################################################
### code chunk number 24: ebnfNonAlternation-repetition (eval = FALSE)
###################################################
## # repetition
## concatenation(
## whitespace(),charParser("{"),whitespace(),
## ebnfDefinition(),
## whitespace(),charParser("}"),
## ebnfConcatenation(),
## action=function (s) 
## if(s[[7]]$value$type=="empty") paste("repetition0N(",s[[4]],")")
## else paste("concatenation(",   paste("repetition0N(",s[[4]],")"),",", s[[7]]$value$value,")",sep=""))


###################################################
### code chunk number 25: ebnfNonAlternation-option (eval = FALSE)
###################################################
## # option
## concatenation(
## whitespace(),charParser("["),whitespace(),
## ebnfDefinition(),
## whitespace(),charParser("]"),
## ebnfConcatenation(),
## action=function (s) 
## if(s[[7]]$value$type=="empty") paste("option(",s[[4]],")")
## else paste("concatenation(",   paste("option(",s[[4]],")"),",", s[[7]]$value$value,")",sep=""))


###################################################
### code chunk number 26: vignette
###################################################
ebnfAlternation <- function() 
concatenation(
ebnfNonAlternation(),

repetition1N(
concatenation(whitespace(),charParser("|"),whitespace(),ebnfNonAlternation(),action=function(s) s[[4]]),action=function(s) s),
action=function(s) paste("alternation(",paste(s[[1]],",",paste(unlist(s[[2]]), collapse=","), sep=""),")",sep=""))
 ebnfNonAlternation <- function() alternation(
# string
concatenation(
string(action=function(s) paste("keyword('",s,"')",sep="")),
ebnfConcatenation(),
action=function (s) 
if(s[[2]]$value$type=="empty") s[[1]]
else paste("concatenation(",s[[1]],",",s[[2]]$value$value,")",sep=""))
                                              ,
# special sequence
concatenation(
ebnfSpecialSequence(),
ebnfConcatenation(),
action=function (s) 
if(s[[2]]$value$type=="empty") s[[1]]
else paste("concatenation(",s[[1]],",",s[[2]]$value$value,")",sep=""))
                                              ,
# rule call
concatenation(
symbolic(charFirst=isLetter,charRest=function(ch) isLetter(ch) || isDigit(ch) || ch == "_",action=function(s) paste(s,"()",sep="")),
ebnfConcatenation(),
action=function (s) 
if(s[[2]]$value$type=="empty") s[[1]]
else paste("concatenation(",s[[1]],",",s[[2]]$value$value,")",sep=""))
                                              ,
# grouping
concatenation(
whitespace(),charParser("("),whitespace(),
ebnfDefinition(),
whitespace(),charParser(")"),
ebnfConcatenation(),
action=function (s) 
if(s[[7]]$value$type=="empty") s[[4]]
else paste("concatenation(",s[[4]],",",s[[7]]$value$value,")",sep=""))
                                              ,
# repetition
concatenation(
whitespace(),charParser("{"),whitespace(),
ebnfDefinition(),
whitespace(),charParser("}"),
ebnfConcatenation(),
action=function (s) 
if(s[[7]]$value$type=="empty") paste("repetition0N(",s[[4]],")")
else paste("concatenation(",   paste("repetition0N(",s[[4]],")"),",", s[[7]]$value$value,")",sep=""))
                                              ,
# option
concatenation(
whitespace(),charParser("["),whitespace(),
ebnfDefinition(),
whitespace(),charParser("]"),
ebnfConcatenation(),
action=function (s) 
if(s[[7]]$value$type=="empty") paste("option(",s[[4]],")")
else paste("concatenation(",   paste("option(",s[[4]],")"),",", s[[7]]$value$value,")",sep=""))
                                              , action=function(s) s)
ebnfDefinition <- function() alternation(
# several alternatives
ebnfAlternation(),
# No alternatives
ebnfNonAlternation(),
action=function(s) s)


###################################################
### code chunk number 27: vignette
###################################################
ebnfSpecialSequence <- function()
concatenation(whitespace(),charParser("?"),whitespace(),
alternation(
keyword("whitespace"      ,action=function(s) s), 
keyword("symbolic"        ,action=function(s) s), 
keyword("string"          ,action=function(s) s), 
keyword("numberInteger"   ,action=function(s) s),
keyword("numberScientific",action=function(s) s),
action=function(s) paste(s,"()",sep="")),
whitespace(),charParser("?"),
action=function(s) s[[4]])


###################################################
### code chunk number 28: vignette
###################################################

  errorFun <- function(strmPosition,h=NULL,type="") { 
    if ( is.null(h) || type != "concatenation" )
      print(paste("Error from line:",strmPosition$line,
      " Caracter:",strmPosition$linePos," Stream Pos:", strmPosition$streamPos, "Type:",type))
    else errorFun(h$pos,h$h,h$type)
      
    return(list(type=type,pos=strmPosition,h=h))
  } 


###################################################
### code chunk number 29: vignette
###################################################
stream  <- streamParserFromString('program = \'PROGRAM\' ;')
cstream <- ebnfRule()(stream)
print(cstream[c("status","node")])	   


###################################################
### code chunk number 30: vignette
###################################################
stream <- streamParserFromString('program =  \'PROGRAM\' , white_space , identifier , white_space ;')
cstream <- ebnfRule()(stream)
print(cstream[c("status","node")])	   


###################################################
### code chunk number 31: vignette
###################################################
stream <- streamParserFromString(
'
program = \'PROGRAM\' , white_space , identifier , white_space ,
          \'BEGIN\'   , white_space ,
           { assignment , ";" , white_space } ,
           \'END.\' ;
')
cstream <- ebnfRule()(stream)
print(cstream[c("status","node")])	   


###################################################
### code chunk number 32: vignette
###################################################
stream <- streamParserFromString(
'identifier = alphabetic_character , { alphabetic_character | digit } ;')
cstream <- ebnfRule()(stream)
print(cstream[c("status","node")])	   


###################################################
### code chunk number 33: vignette
###################################################
stream <- streamParserFromString('white_space = ? whitespace ? ;')
cstream <- ebnfRule()(stream)
print(cstream[c("status","node")])	   


###################################################
### code chunk number 34: vignette
###################################################
stream <- streamParserFromString(
'
program = \'PROGRAM\' , white_space , identifier , white_space ,
          \'BEGIN\'   , white_space ,
           { assignment , ";" , white_space } ,
           \'END.\' ;
identifier = alphabetic_character , { alphabetic_character | digit } ;
number = [ "-" ] , digit , { digit } ;
assignment = identifier , ":=" , ( number | identifier | string_ ) ;
alphabetic_character = "A" | "B" | "C" | "D" | "E" | "F" | "G"
                     | "H" | "I" | "J" | "K" | "L" | "M" | "N"
                     | "O" | "P" | "Q" | "R" | "S" | "T" | "U"
                     | "V" | "W" | "X" | "Y" | "Z" ;
digit = "0" | "1" | "2" | "3" | "4" | "5" | "6" | "7" | "8" | "9" ;
white_space = ? whitespace ? ;
string_ = ? string ?;
')

cstream <- gramatica()(stream)
print(cstream[c("status")])


###################################################
### code chunk number 35: vignette
###################################################
print(cstream[[c("node")]])
eval(parse(text=cstream[[c("node")]]))


###################################################
### code chunk number 36: vignette
###################################################
identifier()(streamParserFromString("DEMO1"))$status

identifier()(streamParserFromString("A0"))$status

keyword(':=')(streamParserFromString(":="))$status

number()(streamParserFromString("3"))$status


###################################################
### code chunk number 37: vignette
###################################################
stream <- streamParserFromString(
'PROGRAM DEMO1
BEGIN
  A0:=3;
  B:=45;
  H:=-100023;
  C:=A;
  D123:=B34A;
  BABOON:=GIRAFFE;
  TEXT:="Hello world!";
END.')


cstream <- program()(stream)
if ( cstream$status=="fail" ) errorFun(cstream$node$pos,cstream$node$h,cstream$node$type) else print(cstream[c("status")])	   



