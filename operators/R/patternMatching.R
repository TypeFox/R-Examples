
### pattern matching
`%~%` <- 
function(x, rx){
  .regexpr <- regexpr %but% getOption("operators.regexpr")
  .regexpr(rx, x) > 0 
}

# inverse of %~%
`%!~%` <- function(x,rx) !`%~%`(x,rx)

# regular expression filters
`%~|%` <- function(x,rx){
  x[x %~% rx]
}
`%!~|%` <- function(x,rx){
  x[! x %~% rx]
}

# all versions                   
`%~*%`  <- function(x,rx) all( `%~%` (x,rx) )
`%!~*%` <- function(x,rx) !all( `%~%`(x,rx) )

# any versions
`%~+%`  <- function(x,rx)  any( `%~%`(x, rx) ) 
`%!~+%` <- function(x,rx) !any( `%~%`(x, rx) )
  

