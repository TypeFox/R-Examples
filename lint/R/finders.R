{############################################################################### 
# finders.R
# This file is part of the R package lint.
# 
# Copyright 2012 Andrew Redd
# Date: 6/16/2012
# 
# DESCRIPTION
# ===========
# Specific region finders.
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

#' @rdname finders
#' @include find-utils.R
#' @title Finder Functions
#' 
#'  Finders are helper functions that assist in the definition of style checks.
#'  These function assist by finding the regions that are either included 
#'  or excluded from check.  Each finder must accept the following arguments
#'  and as well as the variadic argument \code{...}.
#'  \enumerate{
#'      \item file
#'      \item lines
#'      \item parse.data
#'  }
#'  Order of arguments is not guaranteed so explicit names use is required, 
#'  since use of named arguments is guaranteed.
#'  In addition other arguments may be added later.
#'  Each finder is expected to return a \link[lint:conversion]{find} 
#'  formated data.frame.
#'  
#'  Custom finders are encouraged, but finders for common classes are 
#'  included in the package.
#'  
#' @param classes with \code{make_class_finder} the \code{\link{getParseData}}
#'                classes to find.
#' @param file the file being examined.
#' @param lines the lines being examined.
#' @param parse.data data from \code{\link{getParseData}}
#' @param ... extra arguments that include \code{file}, and \code{lines}
#' 
#' @export
make_class_finder <- function(classes){
    f <- function(..., parse.data) {
        rows  <- subset(parse.data, parse.data$token %in% classes)
        if(nrow(rows) == 0) return(empty.find)
        return(rows[c('line1', 'col1', 'line2', 'col2')])
    }
    structure(f, classes=classes)
}


#' @rdname finders
#' @export
find_comment <- make_class_finder(c("COMMENT", "ROXYGEN_COMMENT"))


#' @rdname finders
#' @export
find_basic_comment <- make_class_finder("COMMENT")


#' @rdname finders
#' @export
find_inside_comment <- function(...,parse.data) {
    df <- find_basic_comment(..., parse.data = parse.data)
    df$col1 <- df$col1 + 1
    subset(df, df$col2 > df$col1 | df$line2 > df$line1)
}

#' @rdname finders
#' @export
find_doc_comment <- make_class_finder(c("ROXYGEN_COMMENT"))


#' @rdname finders
#' @export
find_string <- make_class_finder(c("STR_CONST"))

#' @rdname finders
#' @export
find_symbol <- make_class_finder("SYMBOL")

#' @rdname finders
#' @export
find_number <- make_class_finder("NUM_CONST")


#' @rdname finders
#' @export
find_function_args <- function(..., parse.data) {
    ftokens <- subset(parse.data, parse.data$token == "FUNCTION")
    ddply(ftokens, "id" , .find_function_args1
         , parse.data = parse.data)[names(empty.find)]
}
.find_function_args1 <- function(d, ..., parse.data) {
    p <- d$parent
    function.args <- subset(parse.data, parse.data$parent == d$p & 
      !(parse.data$token %in% c('expr', 'FUNCTION')))
    parse2find(function.args)
}
  
#' @rdname finders
#' @export
find_function_body <- function(..., lines, file
    , parse.data = getParseData(parse(file, keep.source=TRUE))) {
  f.nodes <- subset(parse.data, parse.data$token == "FUNCTION")
  if(!nrow(f.nodes)) return(empty.find)
  body.parents  <- ldply(get_children(f.nodes$parent, parse.data, 1), tail, 1)
  body.contents <- find_children(body.parents, parse.data)
  parse2find(body.contents)
}


#' @rdname finders
#' @export
find_call_args <- function(..., file, parse.data = getParseData(parse(file))) {
  call.nodes <- subset(parse.data, 
    parse.data$token == "SYMBOL_FUNCTION_CALL")
  if(!nrow(call.nodes)) return(empty.find)
  call.args <- 
    llply(call.nodes$id, get_family, parse.data=parse.data, nancestors=2)
  parse2find(call.args)
}
