{############################################################################### 
# find-utils.R
# This file is part of the R package lint.
# 
# Copyright 2012 Andrew Redd
# Date: 6/16/2012
# 
# DESCRIPTION
# ===========
# utilities for creating find strip and extract utilities.
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



strip <- function(lines, replace.data, replace.with = ''){
#'  Strip a region from the text
#'  
#'  The \code{strip} function removes the region defined in \code{replace.data} 
#'  from the \code{lines}
#'  
#'  @param lines the lines with the text.  Results from \code{\link{readLines}}
#'  @param replace.data replace data info.  See also \code{\link{find2replace}}
#'  @param replace.with what to replace with , if there is a need.
#'  
#'  @return
#'  The \code{lines} with the regions defined in replace.data removed.
#'  
#'  @family find-functions
#'  @export
  if(nrow(replace.data) == 0) return(lines)
  replace.data <- mutate(replace.data, string = lines[replace.data$line])
  var.names <- c('string', 'start', 'end')
  new.lines <- maply(replace.data[, var.names], `str_sub<-`
                     , value = replace.with, .expand = F)
  lines[replace.data$line] <- new.lines
  lines
}
extract <- function(lines, replace.data) {
#' @rdname strip
#' @description 
#' The \code{extract} function does the opposite of \code{strip}, it extracts
#' the region(s) that were found, dropping everything else.
#' 
#' @export
  if(nrow(replace.data) == 0) return(lines)
  replace.data <- mutate(replace.data, string = lines[replace.data$line])
  var.names <- c('string', 'start', 'end')
  maply(replace.data[, var.names], `str_sub`, .expand=F)
}
make_stripper <- function(finder, replace.with = ''){
  replace.with.default <- replace.with
  function(
    lines,
    text =  paste(lines, collapse='\n'),
    file = textConnection(text), 
    parse.data = getParseData(parse(file, keep.source=TRUE)),
    replace.with = replace.with.default,
    ...
  ){
    find <- finder(parse.data = parse.data)
    strip(lines, find2replace(find), replace.with=replace.with)  
  }
}
make_extractor <- function(finder) {
  function(
    lines,
    text =  paste(lines, collapse='\n'),
    file = textConnection(text), 
    parse.data = getParseData(parse(file, keep.source=TRUE))
  ) {
    find <- finder(parse.data = parse.data)
    extract(lines, find2replace(find))
  }
}
.find_finder_fun1 <- function(ex) {
    if(is.function(ex)) return(ex)
    stopifnot(is.character(ex))
    fun <- try(get(ex, mode = 'function', inherits = T), silent = TRUE)
    if(!inherits(fun, 'try-error')) return(fun)
    fun <- try(match.fun(ex))
    if(!inherits(fun, 'try-error')) return(fun)
    stop(sprintf("Could not find `%s`", ex))
}
find_finder_fun <- Vectorize(.find_finder_fun1, SIMPLIFY = TRUE)

