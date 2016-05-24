{###############################################################################
# test_find.R
# This file is part of the R package lint.
#
# Copyright 2012 Andrew Redd
# Date: 6/16/2012
#
# DESCRIPTION
# ===========
# This file contains the unit tests for the find, strip, extract functions.
#
# LICENSE
# ========
# lint is free software: you can redistribute it and/or modify it under the
# terms of the GNU General Public License as published by the Free Software
# Foundation, either version 3 of the License, or (at your option) any later
# version.
#
# dostats is distributed in the hope that it will be useful, but WITHOUT ANY
# WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
# FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License along with
# this program. If not, see http://www.gnu.org/licenses/.
#
# Log
# ===
# 1/18/2013 tested to work under R-3.0
#
}###############################################################################
context("find/strip/extract")

test_that("Find comments", {
text={"
abc#123
hello # there
"}
  lines <- readLines(textConnection(text))
  comment.locations <- {data.frame(
    line1 = c(2L,  3L),
    col1  = c(4L,  7L),
    line2 = c(2L,  3L),
    col2  = c(7L, 13L),
    stringsAsFactors = FALSE)}

  pd <- getParseData(parse(text=text, keep.source=TRUE))
  expect_that(
      find_comment(parse.data=pd)[,names(comment.locations)]
    , is_equivalent_to(comment.locations))
})
test_that("Find strings", {
  text <- {
'"i\'m a string"
a <- "string"
b <- "second
string"
c <- \'c\'
this.line(has=\'two\', "strings")
no.string'}
  lines <- readLines(textConnection(text))
  pd <- 
  parse.data <- getParseData(parse(text=text, keep.source=TRUE))

  string.loc <- {data.frame(
    line1 = c( 1L,  2L,  3L,  5L,  6L,  6L),
    col1  = c( 1L,  6L,  6L,  6L, 15L, 22L),
    line2 = c( 1L,  2L,  4L,  5L,  6L,  6L),
    col2  = c(14L, 13L,  7L,  8L, 19L, 30L),
    stringsAsFactors = FALSE)}

  expect_that(
    find_string(parse.data=pd)
    , is_equivalent_to(string.loc))
  # results not yet defined
  # expect_that(
  #     strip_string(lines=lines, parse.data=pd)
  #   , equals(c('""', 'a <- ""', 'b <- ""', '""', 'c <- ""'
  #             , 'this.line(has="", "")', "no.string")))
  # expect_that(
  #     extract_string(lines, parse.data=pd)
  #   , equals(c('""', 'a <- ""', 'b <- ""', '""', 'c <- ""', "no.string")))

})
test_that("Find functions arguments", {
  file <- find_example('check-function_args.R', package='lint')
  lines <- readLines(file)
  pd <- 
  parse.data <- getParseData(parse(file, keep.source=TRUE))
  args.loc <- {as.data.frame(structure(matrix(
            c(c( 1L,  9L,  1L, 10L)
            , c( 2L,  9L,  2L, 10L)
            , c( 3L,  9L,  3L, 10L)
            , c( 4L,  9L,  4L, 10L)
            , c( 6L, 14L,  6L, 15L)
            , c( 9L, 14L,  9L, 19L)
            , c(10L, 14L, 10L, 24L)
            , c(13L, 14L, 15L, 18L)
            , c(20L, 20L, 20L, 22L) )
            , byrow=TRUE, ncol=4)
            , dimnames=list(NULL, names(empty.find))))}

  expect_that(
      find_function_args(parse.data=pd, internal=TRUE)
    , equals(args.loc))
  # results not defined
  # expect_that(strip_function_args(lines, parse.data=pd),
  # , equals(lines.no.args))
})
