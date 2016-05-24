# Copyright 2012 Google Inc. All Rights Reserved.
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

Count <- function(x) {
  stopifnot(inherits(x, "histogram"))
  return (sum(x$counts))
}

ApproxMean <- function(x) {
  stopifnot(inherits(x, "histogram"))
  return (weighted.mean(x$mids, x$counts))
}

# TODO(mstokely): Using only the midpoints and counts throws away the
# information about the breakpoints, which could potentially be used
# to provide more accurate approximations.  An alternative slower
# implementation based on cumsum may be added back here.
# Tim suggests using approx there:
# approx(cumsum(c(0, x$counts)) / sum(x$counts, x$breaks, probs, ...)

# wtd.quantile is from Hmisc package.
ApproxQuantile <- function(x, probs, ...) {
  stopifnot(inherits(x, "histogram"))
  if (length(x$mids) < 100) {
    warning("Quantiles computed for histograms with fewer than 100 buckets ",
            "may be inaccurate. Consider using histograms with more granular ",
            "buckets.")
  }
  stopifnot(inherits(x, "histogram"))
  return (wtd.quantile(x$mids, weights=x$counts, probs=probs, ...))
}
