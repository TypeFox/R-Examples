{############################################################################### 
# lint.R
# This file is part of the R package lint.
# 
# Copyright 2012 Andrew Redd
# Date: 6/16/2012
# 
# DESCRIPTION
# ===========
# primary lint functions.
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

#' @name package-lint
#' @title \code{lint}: R code style checker
#' @docType package
#' @author Andrew Redd
#' 
#' @details
#' \code{lint}
#' 
#' @import plyr
#' @import stringr
#' @importFrom harvestr noattr
#' @import foreach
#' @import dostats
#' @include conversion.R
#' @include family.R
#' @include finders.R
NULL

find_region <- function(region, file, lines, parse.data){
    if (length(region) > 0L) {
        fun.region <- find_finder_fun(region)
        if(is.function(fun.region)){
            fun.region(file=file, lines=lines, parse.data=parse.data)
        } else if(is.list(fun.region) && length(fun.region) == 1) {
            fun.region[[1]](file = file, lines = lines, parse.data = parse.data)
        } else if(is.list(fun.region) && length(fun.region) >= 2) {
            l <- vector('list', length(fun.region))
            for(i in seq_along(fun.region))
            l[[i]] <- fun.region[[i]](file = file, lines = lines
                                      , parse.data = parse.data)
            Reduce(merge_find, l)
        } else stop("malformed region!")
    } else empty.find
}

#' Check a pattern against lines
#' 
#' This is part of lint which checks lines of text for stylistic problems.
#' The pattern provided should be a Perl compatible regular expression.
#' 
#' @param lines character vector of lines of text to check, output from 
#'   \code{\link{readLines}}.
#' @param pattern Perl compatible regular expression.
#' @param ... discarded.
#' @return returns an integer vector of the problem lines if problems were 
#'   found. otherwise returns TRUE if the lines pass the check. 
#'   \code{\link{isTRUE}} could be useful for checking the return value.
#' @return \link[lint:conversion]{find} compatible format.
#' @family lint
#' @export
check_pattern <- function(pattern
  , lines
  , ...) {
  if (length(pattern) > 1) {
    pat <- NULL
    foreach(pat=pattern, .combine=merge_find, .multicombine=TRUE
           , .inorder=FALSE) %do% check_pattern(pat, lines, ...)
  } else {
    problem    <- str_locate(pattern=perl(pattern), string=lines)
    locate2find(problem)
  }
}

#' Find an example file.
#' 
#' only for helping with testthat/devtools testing.
#' @keywords internal
find_example <- function(file, package = NULL){
    # cat(getwd(), '\n')
    rf <- system.file("examples", file, package=package, mustWork=FALSE)
    if (rf != "") return(rf) 
    {
        dcf <- file.path('.', 'DESCRIPTION')
        if (file.exists(dcf) && (read.dcf(dcf)[1,'Package'] == package)) {
            rf <- file.path('.', 'inst', 'examples', file)
                if (rf != "") return(rf)
        }
    }
    if(file.exists(package)) {
        dcf <- file.path(package, 'DESCRIPTION')
        if (file.exists(dcf) && (read.dcf(dcf)$Package == package)) {
            rf <- file.path('.', 'inst', 'examples', file)
                if (rf != "") return(rf)
        }
    }
    rf <- file.path('..', 'examples', file)
    if(file.exists(rf)) 
        return(rf)
    stop(sprintf("could not find the file %s", file))
}



#' Look for an argument.
#' 
#' @param x an object
#' @param default the default value
#' @return If x is neither NULL nor NA then x otherwise the default
#' @export
with_default <- function(x, default) {
  if (all(is.null(x))) return(default)
  if (length(x) > 0 && all(is.na(x))) return(default)
  x
}


#' Dispatch tests to the appropriate handler
#' @param test the test
#' @param file the file to check
#' @param parse.data parse data from \code{\link{getParseData}}
#' @param lines the lines to evaluate, overrides file.
#' @param quiet should the test be quiet, i.e. no messages or warnings?
#' @param warning should messages be upgraded to warnings, ignored if 
#'                \code{quiet=TRUE}.
#' 
#' @description
#' runs a test the the appropriate handler.
#' 
#' @return 
#' returns the results from the test handler, which should be either a TRUE for
#' a passed test or the lines, locations, and/or string violating the rules.
#' @export
dispatch_test <- 
function(test, file
        , parse.data = getParseData(parse(file, keep.source=TRUE))
  , lines = readLines(file), quiet = FALSE
  , warning = with_default(test$warning, FALSE)
) {
    include.spans <- find_region(with_default(test$include.region, character(0))
                          , file=file, lines=lines, parse.data=parse.data)
    
    exclude.region <- 
         with_default(test$exclude.region, c("find_comment", "find_string"))
    exclude.spans <- find_region(exclude.region
                          , file=file, lines=lines, parse.data=parse.data)
    
    if(nrow(include.spans) && nrow(exclude.spans))
        stop("specifying both include and exclude regions is undefined")
    
    use.lines <- with_default(test$use.lines, TRUE)
    if (!use.lines) lines <- paste(lines, '\n', collapse='')
    
    do_message <- if(quiet){
        function(...){}
    } else if(warning) {
        get("warning", mode="function")
    } else {
        get("message", mode="function")
    }

    if (!is.null(test$pattern)) {
        test.result <- do.call(check_pattern, append(test, list(lines=lines)))
    } 
    else if(!is.null(test$f)) {
        new.args <- append(test, list(file=file, lines=lines
                                      , parse.data=parse.data))
        test.result <- do.call(check_functional, new.args)
    } 
    else stop("Ill-formatted check.")
    
    if(isTRUE(test.result) || nrow(test.result) == 0) return(TRUE)
    
    # check results.
    if(nrow(exclude.spans) > 0) {
        test.result <- span_difference(test.result, exclude.spans)
    }    
    if(nrow(include.spans) > 0) {
        test.result <- span_intersect(test.result, include.spans)
    }
    if(nrow(test.result)){
        test.message <- with_default(test$message, test$pattern)
        str <- sprintf("Lint: %s: found on lines %s", test.message, 
                       format_problem_lines(test.result$line1))
        do_message(str)
        return(invisible(test.result))
    } else return(TRUE)
}

format_problem_lines <- function(lines, max.to.show = 5) {
    if(length(lines) <= max.to.show) 
        return(paste(lines, collapse = ', '))
    paste(c( head(lines, max.to.show)
           , sprintf("+%d more.", length(lines) - max.to.show))
         , collapse = ', ')
}

   
#' Check for stylistic errors.
#' @param file a vector of file paths.
#' @param style The list of styles tests to use to check.
#' @param text text to check
#' @param recurse recurse into directory and sub-directories.
#' 
#' Check source documents for stylistic errors.  The test are given as a list
#' in  \code{tests}.  If a directory is given all *.R files in that directory 
#' and sub-directories are checked.  If a file other than a .R or .r file 
#' is desired it must be given explicitly as the \code{file} argument.
#' 
#' 
#' @family lint
#' @export
lint <- function(file = '.', style = lint.style, recurse = TRUE, text = NULL) {
    stopifnot(missing(file) | inherits(file, 'character'))
    if (missing(file) && !is.null(text)) {
        stopifnot(inherits(text, "character"))
        file <- textConnection(text, 'rt', local=TRUE)
        files <- list(file)
        on.exit(close(file))
    } else {
        fi <- file.info(file)
        files <- if(any(fi$isdir)) {
            c( file[!fi$isdir]
             , dir(file[fi$isdir], pattern = ".*\\.[Rr]$"
                  , full.names = TRUE, recursive = recurse))
        } else {
            file
        }
    }

    invisible(llply(files, lint_file, style = style))
}

lint_file <- function(file, style) {
    message("Lint checking: ", ifelse(inherits(file, 'character')
                                     , file, class(file)[[1]]))
    lines <- readLines(file)
    parse.data <- getParseData(parse(text=lines, keep.source=TRUE))
  
    invisible(llply(style, dispatch_test, file = file
        , parse.data = parse.data, lines = lines))
}
