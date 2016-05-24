#do not edit, edit noweb/qmrparser.nw
isWhitespace <- function(ch) switch(ch,' '= , '\t'= , '\f'= , '\r'= , '\n'=TRUE,FALSE)
