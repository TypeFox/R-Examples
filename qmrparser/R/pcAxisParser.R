#do not edit, edit noweb/qmrparser.nw
pcAxisParser <- function(streamParser) {
  errorFun <- function(strmPosition,h=NULL,type="") { 
    if ( is.null(h) || type != "concatenation" )
      print(paste("Error from line:",strmPosition$line," position:",strmPosition$linePos))
    else errorFun(h$pos,h$h,h$type)
      
    return(list(type=type,pos=strmPosition,h=h))
  } 
tlistParse <- 
  concatenation(keyword("TLIST"), 
                whitespace(),
                charParser("("),
                whitespace(),
                symbolic(),
                whitespace(),
                alternation(
                            concatenation(charParser(")"),
                                          repetition1N(
                                                       concatenation(
whitespace(), 
                                                                     
charParser(","),
whitespace(),
string(action = function(s) s),
action = function(s) s[[4]]),
                                                       
action = function(s) list(type="list",value=unlist(s))),
                                          action = function(s) s[[2]] ),
                            concatenation(charParser(","),
                                            whitespace(),
                                            string(action = function(s) s),
                                            whitespace(),
                                            charParser("-"),
                                            whitespace(),
                                            string(action = function(s) s),
                                            whitespace(),
                                            charParser(")"),
                                            action = function(s) list(type="limit",value=s[c(3,7)]) ),
                              action = function(s) s),
                  whitespace(),
                  charParser(";"),
                  action = function(s) list(type="tlist",value=s[c(5,7)]) )
##
  rule  <- concatenation(
#                       left rule
#                       keyword
                           symbolic(action = function(s) s),
#                       optional multi-language
                           option(concatenation(
                                                whitespace(),
                                                charParser("["),
                                                whitespace(),
                                                symbolic(action = function(s) s),
                                                whitespace(),
                                                charParser("]"),
                                                action = function(s) s[4]),
                                  action = function(s) if (!is.null(s$type) && s$type=="empty") NULL else s),
#                       optional variables / values
                           option(concatenation(
                                                whitespace(),
                                                charParser("("),
                                                whitespace(),
                                                string(action = function(s) s),
                                                repetition0N(
                                                             concatenation(
whitespace(),
charParser(","),
whitespace(),
string(action = function(s) s),
action = function(s) s[4]), # end concatenation
 action = function(s) if (!is.null(s$type) && s$type=="empty") NULL else s$value),# end repetition0N
                                                whitespace(),
                                                charParser(")"),
                                                action = function(s) if( is.null(s[[5]]) ) s[c(4)] else c(s[c(4)],unlist(s[c(5)]))), # end concatenation
                                  action = function(s) if (!is.null(s$type) && s$type=="empty") NULL else s), # end option
                           whitespace(),
                           charParser("="),
whitespace(),
alternation(
         concatenation(
              string(action = function(s)s),                                                                            repetition1N(concatenation(
                               whitespace(),
                               string(action = function(s) s),
                               action = function(s) s[[2]]),
                             action = function(s) s ),
              whitespace(),
              charParser(";"),
              action = function(s)               list(type="stringstring",value=c(s[c(1)],unlist(s[c(2)]))) ),
concatenation(  string(action = function(s) s)          ,
              repetition0N(
                           concatenation(
                                         whitespace(),
                                         charParser(","),
                                         whitespace(),
                                         string(action = function(s) s),
                                         action = function(s) s[[4]]),
                           action = function(s) if (!is.null(s$type) && s$type=="empty") NULL else s$value),
              whitespace(),
              charParser(";"),
              action = function(s) list(type="liststring",value=if( is.null( s[[2]]  ) ) s[c(1)] else c(s[c(1)],unlist(s[c(2)],use.names = FALSE))) ),
concatenation(  
              alternation(
                          numberScientific(action = function(s) s),
                          string(action = function(s) s),
                          dots  (action = function(s) s),
                          action = function(s) s) ,
              repetition0N(
                            concatenation(
                                         separator(),
                                         alternation(
                                                     numberScientific( action = function(s) s),
                                                     string(action = function(s) s),
                                                     dots  (action = function(s) s),                                                     
                                                     action = function(s) s) ,
                                         action = function(s) s[[2]]),
               action = function(s) if (!is.null(s$type) && s$type=="empty") NULL else s$value), # en repetition0N
              whitespace(),
              alternation(
                          concatenation(charParser(";"),
                                        whitespace(),
                                        option(concatenation( charParser(";"), whitespace())),                                        
                                        option(concatenation( charParser("\032"), whitespace()))
                                       ),
                          concatenation(charParser("\032"),
                                        whitespace()
                                       ),
                          eofMark()),
              action = function(s) list(type="list",value=if( is.null( s[[2]]  ) ) s[c(1)] else c(s[c(1)],unlist(s[c(2)], use.names = FALSE))) ),
tlistParse, 
concatenation(symbolic(action = function(s) s),
              whitespace(),
              charParser(";"),
              action = function(s) list(type="symbol",value=s[[1]])),
              action = function(s) s), ## end alternation
              whitespace(), ## blanks behind ;
                           action = function(s) { rule <- s[c(1,2,3,7)] ; names(rule) <- c("keyword","language","parameters","ruleRight") ; rule }
                           )
cstream <- concatenation(repetition1N(rule,action = function(s) s) ,eofMark(error=errorFun),action = function(s) s[[1]]) (streamParser)
  return(cstream)
}
