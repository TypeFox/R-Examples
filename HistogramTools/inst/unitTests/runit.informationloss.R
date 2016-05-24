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

.checkMetricInvariants <- function(metric) {
  checkTrue(metric >= 0)
  checkTrue(metric <= 1)
}

TestKSDCC <- function() {
  x <- rexp(100)
  h1 <- hist(x, plot=FALSE)
  h2 <- hist(x, breaks=seq(0, round(max(x) + 1), by=0.1), plot=FALSE)

  ksdcc.1 <- KSDCC(h1)
  ksdcc.2 <- KSDCC(h2)

  .checkMetricInvariants(ksdcc.1)
  .checkMetricInvariants(ksdcc.2)
  checkTrue(ksdcc.1 >= ksdcc.2)

  x1.min <- rep(head(h1$breaks, -1), h1$counts)
  x1.max <- rep(tail(h1$breaks, -1), h1$counts)
  checkEquals(unname(ks.test(x1.min, x1.max, exact=F)$statistic), KSDCC(h1))

  x2.min <- rep(head(h2$breaks, -1), h2$counts)
  x2.max <- rep(tail(h2$breaks, -1), h2$counts)
  checkEquals(unname(ks.test(x2.min, x2.max, exact=F)$statistic), KSDCC(h2))
}

TestEMDCC <- function() {
  set.seed(0)
  x <- rexp(100)
  h1 <- hist(x, plot=FALSE)
  h2 <- hist(x, breaks=seq(0, round(max(x) + 1), by=0.1), plot=FALSE)

  emdcc.1 <- EMDCC(h1)
  emdcc.2 <- EMDCC(h2)

  .checkMetricInvariants(emdcc.1)
  .checkMetricInvariants(emdcc.2)
  checkTrue(emdcc.1 >= emdcc.2)

  if (require(emdist)) {
    MinEcdf <- HistToEcdf(h1, f=0)
    MaxEcdf <- HistToEcdf(h1, f=1)

    A1 <- matrix(c(rep(1, length(h1$counts)),
                   h1$mids, MaxEcdf(tail(knots(MinEcdf), -1))), ncol=3)
    A2 <- matrix(c(rep(1, length(h1$counts)),
                   h1$mids, MinEcdf(head(knots(MinEcdf), -1))), ncol=3)
    # emdist seems to use single-precision floating point, thus the need
    # for 2^-23 as the tolerance rather than .Machine$double.eps
    # From http://en.wikipedia.org/wiki/Machine_epsilon
    checkEquals(emd(A1, A2), emdcc.1, tol=2^-23)
  }
}
