{############################################################################### 
# spacing.patterns.R
# This file is part of the R package lint.
# 
# Copyright 2012 Andrew Redd
# Date: 6/16/2012
# 
# DESCRIPTION
# ===========
# utilities for building and checking style checks.
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
#' @title Style checks
#' @name stylechecks
#' @docType data
#' 
#' @aliases 
#'     style_checks
#'     style_tests
#' 
#' @format
#' Each test can be defined in several formats and is very flexible.
#' A test consists of a names list of attributes.
#' \enumerate{
#'   \item \code{pattern} is a pcre compatible \link[base:regex]{regular
#'         expression} that is tested. Can be a character vector of expressions.
#'   \item \code{message} The message to be displayed.
#'   \item \code{include.region} lists regions to restrict the search to.
#'         Can be a character vector specifying the known regions, or a list of 
#'         functions that interpret output from \code{\link{getParseData}}.
#'   \item \code{exclude.region=c('comments', 'string')} lists regions to 
#'         restrict the search to. Operates the sames as \code{include.region}.
#'   \item \code{use.lines=T} should the pattern be evaluated in lines (default)
#'          or as a contiguous block of text.
#'   \item \code{warning=F}
#' }
NULL


#' Test a patten for consistency.
#' 
#' @param check a style check
#' @param ti the test info data
#' @param only.results if true returns results but does not check for
#'                     correspondence.  For debugging.
#' 
#' @return either NULL or throws an error for use with test_that
#' 
test_style <- function(check, ti, only.results = F) {
    p  <- 
    if(is.null(ti$file)) {
        if(is.null(ti$lines)) {
            if(is.null(ti$text))
                stop( paste0("Invalid ti;"
                    , " one of file, lines, or text, must be specified."))
            else {
                ti$lines <- readLines(ti$file)
                parse(text=text, keep.source=TRUE)
            }
        } else {
            if(is.null(ti$text)) {
                ti$text <- paste(ti$lines, collapse='\n')
            } else {
                lines <- readLines(textConnection(ti$text))
                stopifnot(identical(lines, ti$lines))
            }
            parse(text=ti$text, keep.source=TRUE)
        }
    } else {
        if(is.null(ti$lines))
            lines <- readLines(ti$file)
        else
            stopifnot(identical(ti$lines, readLines(ti$file)))
        parse(file=ti$file, keep.source=TRUE)
    }
    pd <- getParseData(p)
    
    results <- suppressMessages(suppressWarnings({
        dispatch_test(check, file = ti$file, lines = ti$lines, parse.data = pd)
    }))
    if(only.results) return(results)
    expect_equivalent(results, ti$results)
    return(invisible(NULL))
}

#' @rdname test_style
#' @param check.name the name of the test as a string.
#' 
#' \code{autotest_style} uses the \code{.testinfo.<<stylename>>} object to 
#' automatically test styles.  The test info object should be a list with 
#' \code{$lines} and \code{$results}. The '\code{$lines}' element is the input 
#' lines and \code{$results} is the find formated data.frame.
#' 
#' @export
autotest_style <- function(check.name, only.results = FALSE) {
    if(!is.character(check.name))
        check.name <- as.character(substitute(c(check.name)))[ - 1]
    
    check <- get(check.name)
    ti <- get(paste0('.testinfo.', check.name))
    if(only.results)
        test_style( check, ti, only.results)
    else 
        test_that(check.name, test_style( check, ti))
}

make_ex_span <- function(line1, col1, line2, col2) {
    data.frame( line1 = line1
              ,  col1 = col1
              , line2 = line2
              ,  col2 = col2)
}

#' @name rr
#' @aliases .rr
#' @title make a results row
#' 
#' A convenience utility for creating the results rows for autotest style test 
#' info
#' @keywords internal
.rr <- function(line1, col1, line2, col2){
    data.frame( line1=line1
              ,  col1=col1
              , line2=line2
              ,  col2=col2  )
}

