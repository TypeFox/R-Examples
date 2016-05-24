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

# Requires David W. Scott's ASH code from ash package on CRAN.
HistToASH <- function(h, m=5, kopt=c(2,2)) {
  # Create an Average Shifted Histogram from a histogram.
  #
  # See ?ash1
  #
  # Args:
  #   h:  An S3 histogram object.
  #   m:  Smooth parameter for ash1 function
  #   kopt: Parameter for ash1 function.
  stopifnot(inherits(h, "histogram"))
  if (!h$equidist) {
    stop("Histogram must be equidist to generate an Average Shifted Histogram")
  }
  # Assuming we have a finely bucketed initial histogram.
  bin1 <- list(nc = h$counts,
               ab = range(h$breaks),
               nskip = 0)
  return(ash1(bin1, m, kopt))
}
