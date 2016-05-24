# Copyright 2014 Google Inc. All rights reserved.
#
# Use of this source code is governed by a BSD-style
# license that can be found in the LICENSE file or at
# http://opensource.org/licenses/BSD-3-Clause

# R does not have colVars or colStdevs. Create quick versions.

colVars <- function(x, na.rm = FALSE) {
  # Compute variance of each column of x
  # Quick version, does not support dims.
  if (is.data.frame(x)) {
    f <- function(v, na.rm = na.rm) {
      if(is.numeric(v) || is.logical(v) || is.complex(v))
        var(v, na.rm = na.rm)
      else NA # Just return NA even if v is a character matrix
    }
    return(unlist(lapply(x, f, na.rm = na.rm)))
  }
  if(!is.matrix(x))
    return(var(x, na.rm = na.rm))

  means <- colMeans(x, na.rm = na.rm)
  n <- length(x)/length(means)
  centeredX <- x - rep(means, each = n)
  colSums(centeredX^2, na.rm = na.rm) /
    pmax(0, IfElse(anyNA(x), n - 1 - colSums(is.na(x)), n - 1))
}

colStdevs <- function(x, ...)
  sqrt(colVars(x, ...))
