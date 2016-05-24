# Copyright 2013 Google Inc. All Rights Reserved.
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.
#
# Author: mstokely@google.com (Murray Stokely)

HistToEcdf <- function(h, method="constant", f=0, inverse=FALSE) {
  # Compute an empirical cumulative distribution function from a histogram.
  #
  # Args:
  #   h:  An S3 histogram object.
  #   method: specifies the interpolation method to be used in call to
  #      approxfun().  Choices are "linear" or "constant".
  #   f: for method="constant" a number between 0 and 1 inclusive,
  #      indicating a compromise between left- and right-continuous
  #      step functions.  See ?approxfun
  #   inverse: If TRUE, return the inverse function.
  #
  # Returns:
  #   An empirical cumulative distribution function (see ?ecdf)

  n <- sum(h$counts)
  # We don't want to use h$mids here, because we want at least
  # to get the correct answers at the breakpoints.
  x.vals <- h$breaks
  y.vals <- c(0, cumsum(h$counts) / n)
  if (inverse) {
    vals.tmp <- x.vals
    x.vals <- y.vals
    y.vals <- vals.tmp
  }
  rval <- approxfun(x.vals, y.vals,
                    method = method, yleft = 0, yright = 1, f = f,
                    ties = "ordered")
  class(rval) <- c("ecdf", "stepfun", class(rval))
  assign("nobs", n, envir = environment(rval))
  attr(rval, "call") <- sys.call()
  rval
}
