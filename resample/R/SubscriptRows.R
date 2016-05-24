# Copyright 2014 Google Inc. All rights reserved.
#
# Use of this source code is governed by a BSD-style
# license that can be found in the LICENSE file or at
# http://opensource.org/licenses/BSD-3-Clause

# This is based on SubsetRows from library(dataframe)
# http://cran.fhcrc.org/web/packages/dataframe

.resampleSubscriptRows <- function(x, i) {
  # Subscript rows of a data frame.
  #
  # Args:
  #   x: a data frame
  #   i: numeric, rows to extract
  #
  # .resampleSubscriptRows(x, i) is equivalent to x[i, , drop = FALSE]
  # except that the result has artificial row names,
  # like those resulting from 'data.frame(..., row.names = NULL)'.
  # This is faster, because it does not check for
  # and eliminate duplicate row names.
  # This is optimized for speed in other ways.
  y <- vector("list", length(x))
  for(j in seq_along(y)) {
    xj <- .subset2(x, j) # x[[j]]
    y[[j]] <- if(length(dim(xj)) != 2L) xj[i] else xj[i, , drop = FALSE]
  }
  names(y) <- names(x)
  # compute nrows, using a variable if there is one, otherwise length((1:n)[i])
  nrows <- IfElse(length(x) == 0, length(seq_len(nrow(x))[i]),
                  length(dim(y[[1]]) == 2), nrow(y[[1]]),
                  length(y[[1]]))
  attr(y, "row.names") <- .set_row_names(nrows)
  class(y) <- "data.frame"
  y
}


if(FALSE)
.resampleSubscriptRowsColumns <- function(x, i, j) {
  # Subscript rows, of some columns, of a data frame.
  # This is used for permuting selected columns.
  # The whole data frame is returned, with other columns unchanged.
  #
  # Args:
  #   x: a data frame
  #   i: numeric, rows to subscript. This must subscript as many rows
  #      as the original x.
  #   j: columns to subscript.

  # Convert j to positive numerical values
  if(is.character(j))
    j <- match(j, names(x))
  if(any(j < 0) || is.logical(j))
    j <- seq_along(x)[j]
  if(anyNA(j))
    stop("Unable to figure out resampleColumns")
  j <- j[j != 0]

  if(!length(j))
    stop("No columns selected using resampleColumns")

  n <- nrow(x)
  y <- x
  for(jj in j) {
    xj <- .subset2(x, jj) # x[[jj]]
    y[[jj]] <- if(length(dim(xj)) != 2L) xj[i] else xj[i, , drop = FALSE]
  }

  # Check that the number of rows is right (using one variable)
  j1 <- j[[1]]
  nObs <- IfElse(length(dim(y[[j1]]) == 2), nrow(y[[j1]]), length(y[[j1]]))
  if(nObs != n)
    stop("When using resampleColumns you must select n rows.")
  y
}
# TODO: finish this (attributes, class) and use it in MakeFunction.R
# for the resampleColumns case.
