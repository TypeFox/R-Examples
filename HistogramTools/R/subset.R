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

# Subset histogram should return the parts of the histogram between the
# specified breakpoints minbreak and maxbreak.
SubsetHistogram <- function(x, minbreak=NULL, maxbreak=NULL) {
  stopifnot(inherits(x, "histogram"))
  if (is.null(minbreak) && is.null(maxbreak)) {
    warning("No new breakpoints specified, returning original histogram.")
    return(x)
  }
  if (!is.null(minbreak)) {
    stopifnot(is.numeric(minbreak), length(minbreak) == 1)
    stopifnot(minbreak %in% x$breaks)
    if (minbreak == min(x$breaks)) {
      return(x)
    }
    # How many bins to cut from left side of histogram?
    num.to.cut <- sum(x$breaks < minbreak)
    x$breaks <- tail(x$breaks, -num.to.cut)
    x$counts <- tail(x$counts, -num.to.cut)
  }
  if (!is.null(maxbreak)) {
    stopifnot(is.numeric(maxbreak), length(maxbreak) == 1)
    stopifnot(maxbreak %in% x$breaks)
    if (maxbreak == max(x$breaks)) {
      return(x)
    }
    # How many bins to cut from right side of histogram?
    num.to.cut <- length(which(x$breaks > maxbreak))
    x$breaks <- head(x$breaks, -num.to.cut)
    x$counts <- head(x$counts, -num.to.cut)
  }
  # Return a new histogram to update density, mids fields.
  return(.BuildHistogram(x$breaks, x$counts, x$xname))
}
