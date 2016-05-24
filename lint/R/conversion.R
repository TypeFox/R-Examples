{###############################################################################
# Conversion.R
# This file is part of the R package lint.
#
# Copyright 2012 Andrew Redd
# Date: 6/16/2012
#
# DESCRIPTION
# ===========
# functions for conversion between the difference formats
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

#' @name conversion
#' @title Lint internal data structures
#' @keywords utils, internal
#'
#' @section Introduction:
#'  lint makes use of several functions from different packages that
#'  store data in various different formats.  These functions provide
#'  utilities for converting between the different formats.
#'
#'  The formats are:
#'  \enumerate{
#'      \item parse   - from the \code{\link{getParseData}} function.
#'                      In parse data each element of an expression
#'                      has it's own row.
#'      \item find    - similar to parse but gives a row for each region or
#'                      expression of interest.
#'      \item replace - for use with \code{\link[stringr:str_sub]{stringr}}.
#'                      Uses a column structure with start and end, organized
#'                      into a matrix with a row for each line.
#'      \item locate  - results from \code{\link{str_locate}} from stringr.
#'                      same as replace for most purposes but does not include
#'                      a string.
#'  }
#'
#'  @section parse data structure:
#'   Parse data structure originates from the \code{\link{getParseData}}
#'   function,which returns an objects with the attribute '\code{data}'.
#'   Parse formatteddata contains a row for every token, string, and expression.
#'   The data frame describes a tree structure with each row a node.
#'   Each node has a parent unless it is a root node i.e. parent==0.
#'   It has the following columns.
#'   \enumerate{
#'       \item \code{line1} starting line of the expression.
#'       \item \code{col1} starting column.
#'       \item \code{line2} ending line of the expression.
#'       \item \code{col2} ending column.
#'       \item \code{token} the token class number.
#'       \item \code{id} the unique id of the expression
#'       \item \code{parent} the parent of the expression, 0 if none.
#'       \item \code{top_level} top_level, which top level expression is the
#'                              expression associated with
#'       \item \code{token} class name of the token.
#'       \item \code{terminal} is this a terminal node? i.e. has no child nodes.
#'       \item \code{text}  the actual text of the expression.
#'   }
#'
#'   The parse data is formatted with C based indexing.  E.g. the first
#'   two elements would be listed as \code{col1=0, col2=2}.  The line number
#'   however is 1 based so the first line is 1, there is no zero line.
#'
#'
#'  @section Find data structure:
#'   For the purposes of the data the find data consists of a single row
#'   for each section/region with the first 4 columns of parse.data;
#'   the columns \code{line1}, \code{col1}, \code{line2}, and
#'   \code{col2}, marking the beginning and end of a section.
#'   This is a condensation of the parse data which would have the same columns
#'   as well as additional columns, and a row for each expression in the region.
#'
#'   Find formatted data is defined to be R or 1 based arrays and inclusive.
#'   the first two elements would be \code{col1=1, col2=2}.
#'
#'   Although both \code{col} elements are retained in
#'   conversion functions, at this time only \code{col} columns are used
#'   internally.
#'
#'
#'  @section Replace data structure:
#'   The data structure for replace data is defined as a data frame with
#'   columns suitable for use ase arguments to str_sub.  That is it has columns
#'   \enumerate{
#'     \item \code{start}
#'     \item \code{end}
#'     \item and either \code{string} or \code{line}
#'   }
#'   where \code{string} would be preferred but line to match up with line data.
#'   \code{find2replace} uses the line, since the string is not available in the
#'   find data.
#'
#'   Replace data formatted data is also R/1 inclusive based arrays.
#'
#'
#' @section Locate data structure:
#'  locate data is defined as the matrix that comes from
#'  \code{\link{str_locate}}.
#'  It has columns
#'  \enumerate{
#'      \item \code{start}
#'      \item \code{end}
#'  }
#'  and has a row for every line.
NULL

#' @rdname conversion
#' @export
empty.find <- {data.frame(
      'line1' = integer(0L)
    , 'col1'  = integer(0L)
    , 'line2' = integer(0L)
    , 'col2'  = integer(0L) )}

#' @rdname conversion
#' @details parse2find
#'  Deprecated.  Due to the changes in R version 3.0 this function is no longer
#'  necessary
#'
#'  Expects either a parse formatted data.frame or a list of data.frames.
#'  each data.frame is a contiguous region that is collapsed into a single
#'  find formatted data.frame, one row for each region.
#'
#' @export
parse2find <- function(parse.data) {
  col1 <- NULL
  if (!inherits(parse.data, 'data.frame') && inherits(parse.data, 'list')) {
    return(ldply(parse.data, parse2find)[names(empty.find)])
  }
  if(nrow(parse.data) == 0) return(empty.find)
  names1 <- c('line1', 'col1')
  names2 <- c('line2', 'col2')
  pd1 <- parse.data[do.call(order, parse.data[names1]), ]
  pd2 <- parse.data[do.call(order, parse.data[names2]), ]
  mutate(data.frame(
      head(pd1, 1L)[names1]
    , tail(pd2, 1L)[names2]
  ))
}

#' @rdname conversion
#' @export
find2replace <- function(find.data) {
    mdply(find.data, .find2replace1)[, - seq_len(ncol(find.data))]
}
.find2replace1 <- function(line1, col1, line2, col2){
    if (line1 == line2) {
      data.frame(line = line1, start = col1, end = col2)
    } else {
      nlines <- line2 - line1
      data.frame(
        line  = c(line1:line2),
        start = c(col1, rep(1L, nlines)),
        end   = c(rep( - 1L, nlines), col2))
    }
}

#' @rdname conversion
#' @export
locate2find <- function(loc) {
    if(!inherits(loc, 'data.frame'))
        loc <- as.data.frame(loc)
    if(all(is.na(loc$start)))return(empty.find)
    mutate(loc
      , line1 = seq_along(start)
      , line2 = seq_along(start)
      ,  col1 = start
      ,  col2 = end)[!is.na(loc$start), names(empty.find)]
}

do_results_overlap_1 <- function(x, y, strict.contains) {
    #' @note assumes that x and y are 1 row each.
    if (x$line2 < y$line1) return(FALSE)
    if (x$line1 > y$line2) return(FALSE)
    x.start <- x$line1
    x.end   <- x$line2
    y.start <- y$line1
    y.end   <- y$line2
    max.col <- max(x$col1, x$col2, y$col1, y$col2)
    if (max.col > 0) {
        x.start <- x.start + x$col1 / max.col
        x.end   <- x.end   + x$col2 / max.col
        y.start <- y.start + y$col1 / max.col
        y.end   <- y.end   + y$col2 / max.col
    }
    if(strict.contains) {
        if (x.start >= y.start && x.end <= y.end) return(TRUE)
    } else {
        if (x.start <= y.start && y.start <= x.end) return(TRUE)
        if (x.start <= y.end   && y.end   <= x.end) return(TRUE)
        if (y.start <= x.start && x.start <= y.end) return(TRUE)
        if (y.start <= x.end   && x.end   <= y.end) return(TRUE)
    }
    return(FALSE)
}
do_results_overlap <- function(x, y = x, strict.contains = FALSE) {
    #' @note assumes x and y are find results formatted data frames.
    force(x)
    force(y)
    y <- mlply(y, data.frame)
    x <- mlply(x, data.frame)
    z <- matrix(NA, length(x), length(y))
    for(i in seq_along(x)) for(j in seq_along(y))
        z[i, j] <- do_results_overlap_1(x[[i]], y[[j]]
                                        , strict.contains = strict.contains)
    return(z)
    #' @return logical matrix of dimension \code{nrow(x)} by \nrow{nrow(y)}.
}

merge_find <- function(...){
  # merge multiple find results
  # @param ... find results.
  find.results <- list(...)
  if(length(find.results) == 0)return(empty.find)
  else if(length(find.results) == 1){
    if(valid_find(find.results[[1]], T, T))
      return(find.results[[1]])
    else stop("argument is not a strict find result.")
  }
  else if(length(find.results) > 2){
    return(Reduce(merge_find, find.results))
  }
  else {
    x <- find.results[[1]]
    y <- find.results[[2]]
    if(nrow(y) == 0) return(x)
    if(nrow(x) == 0) return(y)
    keep <- names(empty.find)
    overlaps <- data.frame(which(do_results_overlap(x,y), arr.ind = T))
    names(overlaps) <- c('x.idx', 'y.idx')
    if(nrow(overlaps) >= 1) {
        merged <- mdply(overlaps, .merge_by_idx1, x = x, y = y)
        new.finds <- rbind(merged[keep],
        x[ - overlaps$x.idx, keep],
        y[ - overlaps$y.idx, keep])
        new.finds[do.call(order, new.finds), ]
    } else {
        combined <- rbind(x[keep], y[keep])
        combined[do.call(order, combined),]
    }
    #' @results a single \code{\link{data.frame}} with find results where
    #'  overlaps were merged
  }
}
.merge_by_idx1 <- function(x.idx, y.idx, x, y){
    x.row <- x[x.idx, ]
    y.row <- y[y.idx, ]
    cbind( rename(min_find_2( x.row[c('line1', 'col1')]
                            , y.row[c('line1', 'col1')]), .names1)
         , rename(max_find_2( x.row[c('line2', 'col2')]
                            , y.row[c('line2', 'col2')]), .names2) )
}
.names1 <- c(line='line1', col='col1')
.names2 <- c(line='line2', col='col2')
min_find <- function( a.line, a.col
                    , b.line, b.col)
{
    if(a.line == b.line){
        data.frame( line = a.line
                  , col  = min(a.col, b.col) )
    } else
    if(a.line < b.line){
        data.frame( line = a.line
                  , col  = a.col )
    } else
    if(a.line > b.line) {
        data.frame( line = b.line
                  , col  = b.col  )
    } else stop("Seriously?! how did you get here?")
}
min_find_2 <- function(x, y)
    min_find(x[['line1']], x[['col1']], y[['line1']], y[['col1']])
max_find <- function( a.line, a.col
        , b.line, b.col) {
    if(a.line == b.line){
        data.frame( line = a.line
                  , col  = max(a.col, b.col) )
    } else
    if(a.line > b.line){
        data.frame( line = a.line
                  , col  = a.col  )
    } else
    if(a.line < b.line) {
        data.frame( line = b.line
                  , col  = b.col  )
    } else stop("Seriously!? how did you get here?")
}
max_find_2 <- function(x, y)
    max_find(x[['line2']], x[['col2']], y[['line2']], y[['col2']])




valid_find <- function(x, strict = FALSE, extended = TRUE){(
    is(x,'data.frame')
 && if(strict) {
        identical(names(empty.find), names(x))
    } else {
        all(names(empty.find) %in% names(x))
    }
 && if(extended) {
        !any(do_results_overlap(x) & !diag(T, nrow(x), nrow(x)))
    } else {
        TRUE
    }
 && !any(x$line1 == 0) && !any(x$col1 == 0)
 && !any(x$line1 < x$line2)
)}

span_intersect <- function(x, y){
    if(nrow(x) == 0 || nrow(y) == 0) return(empty.find)
    overlaps <- data.frame(which(do_results_overlap(x,y), arr.ind=T))
    if(nrow(overlaps) == 0) return(empty.find)
    names(overlaps) <- c('x.idx', 'y.idx')
    mdply(overlaps, .span_intersect1, x=x, y=y)[names(empty.find)]
}
.span_intersect1 <- function(x.idx, y.idx, x, y){
    x.row <- x[x.idx, ]
    y.row <- y[y.idx, ]

    line1 <- max(x.row$line1, y.row$line1)
    col1  <- if(x.row$line1 == y.row$line1) max(x.row$col1, y.row$col1)
             else if(x.row$line1 < y.row$line1) y.row$col1
             else x.row$col1
    line2 <- min(x.row$line2, y.row$line2)
    col2  <- if(x.row$line2 == y.row$line2) min(x.row$col2, y.row$col2)
             else if(x.row$line2 > y.row$line2) y.row$col2
             else x.row$col2

    data.frame( line1 = line1
              ,  col1 =  col1
              , line2 = line2
              ,  col2 =  col2)
}

span_difference <- function(x, y, strict = FALSE) {
#' @param x find formatted data.frame
#' @param y find formatted data.frame
#' @param strict should a strict difference be taken or only
#'               those of x not completely contained in y (default).
    if(nrow(x) == 0 || nrow(y) == 0) return(x)
    overlaps <- do_results_overlap(x,y, TRUE)
    xo <- apply(overlaps, 1, any)
    if(!strict) {
        return(x[!xo,])
    } else stop('I have not done strict span difference yet, not sure I will.')
}

