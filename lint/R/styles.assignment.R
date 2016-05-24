{############################################################################### 
# styles.assignment.R
# This file is part of the R package lint.
# 
# Copyright 2012 Andrew Redd
# Date: 6/16/2012
# 
# DESCRIPTION
# ===========
# predefined assignment styles.
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

#' @title Assignment Style
#' @name assignment-styles
#' @docType data 
#' @format style tests for assignment operators.
#'  
#' @include finders.R
#' @exportPattern styles\\.assignment\\..*
NULL

#' @rdname assignment-styles
#' @aliases styles.assignment.noeq
styles.assignment.noeq <- {list(
    f = make_class_finder('EQ_ASSIGN')
  , message = "Equal sign assignments"
  , exclude.region = .no.exclude
)}
.testinfo.styles.assignment.noeq <- {list(
    lines = c('a=1', 'a <- 1', '1 -> 1', 'f(a=1)', 'f((a=1))')
  , results = data.frame( line1 = c(1, 5)
                        ,  col1 = c(2, 5)
                        , line2 = c(1, 5)
                        ,  col2 = c(2, 5) )
)}
 
#' @rdname assignment-styles
styles.assignment.norightassign <- {list(
    f = make_class_finder('RIGHT_ASSIGN')
  , message = "Right assignment not allowed"
  , exclude.region = .no.exclude
)}
.testinfo.styles.assignment.norightassign <- {list(
    lines = c('a=1', 'a <- 1', '1 -> 1', 'f(a=1)', 'f((a->1))')
  , results = data.frame( line1 = c(3, 5)
                        ,  col1 = c(3, 5)
                        , line2 = c(3, 5)
                        ,  col2 = c(4, 6) )
)}
 
#' @rdname assignment-styles
styles.assignment.notinfcall <- {list(
    f = make_class_finder(c('LEFT_ASSIGN', 'EQ_ASSIGN', 'RIGHT_ASSIGN'))
  , message = "Assignments are not allowed to be nested in function calls."
  , include.region = 'find_call_args'
  , exclude.region = .no.exclude
)}
.testinfo.styles.assignment.notinfcall <- {list(
    lines = c( 'a=1', 'a <- 1', '1 -> a', 'f(a=1)'  # Pass
             , 'f((a=1))', 'f(a <- 1)', 'f(1 -> a)')   # Fail
  , results = data.frame( line1 = c(5, 6, 7)
                        ,  col1 = c(5, 5, 5)
                        , line2 = c(5, 6, 7)
                        ,  col2 = c(5, 6, 6))
)}


