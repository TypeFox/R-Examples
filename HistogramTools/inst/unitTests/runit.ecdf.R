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

TestHistToEcdfInverse <- function() {
  set.seed(0)
  x <- rexp(100)
  h <- hist(x, plot=FALSE)
  f <- HistToEcdf(h)
  f.inv <- HistToEcdf(h, inverse=TRUE)
  checkIdentical(knots(f.inv),
                 sapply(knots(f.inv), function(x) f(f.inv(x))))
}

TestHistToEcdf <- function() {
  set.seed(0)
  x <- rexp(100)
  h <- hist(x, plot=FALSE)
  f <- HistToEcdf(h)
  checkTrue(min(knots(f)) < min(range(x)))
  checkTrue(max(knots(f)) > max(range(x)))

  checkEquals(f(h$breaks[3]), (cumsum(h$counts) / sum(h$counts))[2])
}
