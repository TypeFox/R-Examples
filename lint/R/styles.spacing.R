{############################################################################### 
# styles.spacing.R
# This file is part of the R package lint.
# 
# Copyright 2012 Andrew Redd
# Date: 6/16/2012
# 
# DESCRIPTION
# ===========
# predefined spacing patterns.
# 
# LICENSE
# ========
# lint is free software: you can redistribute it and/or modify it under the
# terms of the GNU General Public License as published by the Free Software 
# Foundation, either version 3 of the License, or (at your option) any later 
# version.
# 
# lint is distributed in the hope that it will be useful, but WITHOUT ANY 
# WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS 
# FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
# 
# You should have received a copy of the GNU General Public License along with 
# this program. If not, see http://www.gnu.org/licenses/.
# 
}###############################################################################

#' @include pattern-utils.R
#' @exportPattern ^spacing\\..*$
NULL


#' @rdname stylechecks
#' @export
spacing.linelength.80 <- {list(pattern = "^.{80}\\s*[^\\s]"
  , message = "Line width exceeds 80 characters"
  , use.lines = TRUE
  , exclude.region = .no.exclude
)}
.testinfo.spacing.linelength.80 <- {list(
    lines = c('123'
      , paste(rep.int('#', 80), collapse='')
      , paste(rep.int('#', 81), collapse=''))
  , results = data.frame(line1=3, col1=1, line2=3, col2=81)
)}

#' @rdname stylechecks
#' @export
spacing.linelength.100 <- {list(pattern = "^.{100}\\s*[^\\s]"
  , message = "Line width exceeds 100 characters"
  , use.lines = TRUE
  , exclude.region = .no.exclude
  , warning = TRUE
)}
.testinfo.spacing.linelength.100 <- {list(
    lines = c('123'
      , paste(rep.int('#', 100), collapse='')
      , paste(rep.int('#', 101), collapse=''))
  , results = data.frame(line1=3, col1=1, line2=3, col2=101)
)}

#' @rdname stylechecks
#' @export
spacing.indentation.notabs <- list(pattern ="^\\t"
  , message = "tabs not allowed for indentation"
  , exclude.region = .no.exclude
)
.testinfo.spacing.indentation.notabs <- {list(
    lines = c('    "hello world"'    #  Good
            , '\t"hi there"'         #  Bad
            , '"don\'t\t catch me"'  #  Good, inside string
            , '#don\'t\t catch me'   #  OK, not at beginning
            , 'IM <-\twrong()')      #  OK, not at beginning
  , results = data.frame(
        line1 = as.integer(c(2))
      , col1  = as.integer(c(1))
      , line2 = as.integer(c(2))
      , col2  = as.integer(c(1)) )
)}

#' @rdname stylechecks
#' @export
spacing.notabs <- list(pattern = "\\t"
  , message = "tabs not ever allowed"
  , exclude.region = "find_string"
)
.testinfo.spacing.notabs <- {list(
    lines = c('    "hello world"'    #  Good
            , '\t"hi there"'         #  Bad
            , '"don\'t\t catch me"'  #  OK, inside string(excluded)
            , "#don't\t catch me"    #  Bad, inside comment, not excluded
            , 'IM <-\twrong()')      #  Bad
  , results = data.frame(
        line1 = as.integer(c(2, 4, 5))
      , col1  = as.integer(c(1, 7, 6))
      , line2 = as.integer(c(2, 4, 5))
      , col2  = as.integer(c(1, 7, 6)) )
)}

#' @rdname stylechecks
#' @export
spacing.indentation.evenindent <- list(pattern = "^(  )*( )\\S"
  , message = "indentation should be by two spaces."
  , exclude.region = c("find_function_args", "find_call_args")
)

#' @rdname stylechecks
#' @export
spacing.spaceaftercomma <- list(pattern = ",[^\\s]"
  , message =  "no space after comma")

#' @rdname stylechecks
#' @export
spacing.spacearoundinfix <- {list(
    pattern = c(paste0(no.preceeding.space.rx, no.preceeding.percent
                      , '(', infix.noeq, ')')
              , paste0( '(', infix.noeq, ')'
                      , no.trailing.percent, no.trailing.space.rx))
  , message = "needs space around infix operators"
  , exclude.region = c("find_comment", "find_string"
                      , "find_symbol", "find_number")
)}
.testinfo.spacing.spacearoundinfix <- {list(
    lines = {c( '1 + 1'                #   1 #
             , '1+2'                   #     #
             , '1+ 2'                  #     #
             , '1 +2'                  #     #
             , '1-2'                   #   5 #
             , '1*2'                   #     #
             , '1/2'                   #     #
             , '1%/%2'                 #     #
             , '1^2'                   #     #
             , '1**2'                  #  10 #
             , '1%*%2'                 #     #
             , '1%o%2'                 #     #
             , '1%in%2'                #     #
             , 'T&F'                   #     #
             , 'T&&F'                  #  15 #
             , 'T|F'                   #     #
             , 'T||F'                  #     #
             , '1==2'                  #     #
             , '1>=2'                  #     #
             , '1<=2'                  #  20 # 
             , '1> 2'                  #     #
             , '1< 2'                  #     #
             , '1!=2'                  #     #
             , 'if(a==b)return(TRUE)'  #     #
             , 'r1 <- x %*% y'         #  25 #  OK
             , 'r2 <- x %/% y'         #     #  OK
             , 'r3 <- x %o% y'         #     #  OK
             , 'r4 <- x %in% y'        #     #  OK
             , 'r5 <- x %myopp% y'     #     #  OK
             , 'r6 <- x %% y'          #  30 #  OK
             , 'a <- 1e-3'             #     #  OK
             , 'b <- 2E+4'             #     #  OK
             , '1 +'                   #     #  OK 
             , '1 + 1'                 #     #  OK 
             )}
  , results = {rbind.fill(
#     line1, col1, line2, col2
  .rr(    2,    2,     2,    2)
, .rr(    3,    2,     3,    2)
, .rr(    4,    3,     4,    3)
, .rr(    5,    2,     5,    2)
, .rr(    6,    2,     6,    2)
, .rr(    7,    2,     7,    2)
, .rr(    8,    2,     8,    4)
, .rr(    9,    2,     9,    2)
, .rr(   10,    2,    10,    3)
, .rr(   11,    2,    11,    4)
, .rr(   12,    2,    12,    4)
, .rr(   13,    2,    13,    5)
, .rr(   14,    2,    14,    2)
, .rr(   15,    2,    15,    3)
, .rr(   16,    2,    16,    2)
, .rr(   17,    2,    17,    3)
, .rr(   18,    2,    18,    3)
, .rr(   19,    2,    19,    3)
, .rr(   20,    2,    20,    3)
, .rr(   21,    2,    21,    2)
, .rr(   22,    2,    22,    2)
, .rr(   23,    2,    23,    3)
, .rr(   24,    5,    24,    6)
)}
)}

#' @rdname stylechecks
#' @export
spacing.spacearoundequals <- {list(
    pattern = c(paste0(no.preceeding.space.rx, '(?<![=!<>])(=)(?![=])')
              , paste0('(?<![=!<>])(=)(?![=])', no.trailing.space.rx))
  , message = "needs space around `=`"
  , exclude.region = c("find_call_args", "find_comment", "find_string")
)}
.testinfo.spacing.spacearoundequals <- {list(
    lines = c( 'a=1'                     #  Bad
             , 'f(a=1)'                  #  ok
             , 'function(){'             #  
             , '    a=1'                 #  bad
             , '}'                       #  
             , 'print(paste("x =", x))'  #  ok
             , 'print(paste("x =",'      #  ok
             ,  'x))'                    #  
             )
  , results = data.frame( line1 = c(1, 4)
                        ,  col1 = c(2, 6)
                        , line2 = c(1, 4)
                        ,  col2 = c(2, 6) )
)}



#' @rdname stylechecks
#' @export
spacing.twobeforecomments <- {list(
    pattern = perl("(?<=[^\\s\\{\\}#])(#)|(?<=[^\\s]\\s)(#)")
  , exclude.region = c("find_string", 'find_inside_comment')
  , message = "needs two spaces spacing before inline comments")
}
.testinfo.spacing.twobeforecomments <- {list(
    lines = {c(  '{#'         #  OK
               , '}#'         #  OK
               , '# c'        #  OK, start of line
               , '1#c'        #  BAD
               , '1 #c'       #  BAD
               , '1 <- "#c"'  #  OK, in string
               , '1  # c #'   #  OK, inside other comment
               )}
  , results = {data.frame(  line1 = c(4, 5)
                          ,  col1 = c(2, 3)
                          , line2 = c(4, 5)
                          ,  col2 = c(2, 3) )}
)}
