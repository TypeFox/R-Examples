# Copyright 2014 Google Inc. All rights reserved.
#
# Use of this source code is governed by a BSD-style
# license that can be found in the LICENSE file or at
# http://opensource.org/licenses/BSD-3-Clause

cat0 <- function(...) cat(..., sep = "")

catn <- function(...) {
  # cat, but with a final newline.
  # catn("text") is equivalent to cat("text\n")
  if (any(nzchar(names(match.call())))) {
    # Quick and dirty way to handle file and append arguments;
    # this results in an extra " " before the final "\n".
    cat(..., "\n")
  } else {
    cat(...)
    cat("\n")
  }
}

cat0n <- function(...) {
  # cat(), but with sep = "" and a final newline.
  # cat0n("a", "b") is equivalent to cat("a", "b", "\n", sep = "")
  cat(..., sep = "", "\n")
}
