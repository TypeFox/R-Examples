# Copyright 2012 Google Inc. All Rights Reserved.
# Author: mstokely@google.com (Murray Stokely)

TestApproxMean <- function() {
  x <- hist(c(1, 2, 3, 4), plot=FALSE)
  # x$breaks = 1 2 3 4
  # x$counts = 2 1 1
  checkEquals(ApproxMean(x), 2.25)
}

TestApproxQuantile <- function() {
  x <- hist(c(1, 2, 3, 4), plot=FALSE)
  # The bucketing above means a reasonable approximation for the median
  # would be between 1.5 and 2.5.
  checkTrue(unname(ApproxQuantile(x, .5)) >= 1.5)
  checkTrue(unname(ApproxQuantile(x, .5)) <= 2.5)
  checkTrue(ApproxQuantile(x, .5) < ApproxQuantile(x, .51))
  # The midpoint of the last bucket is our approximation of 100%tile.
  checkEquals(unname(ApproxQuantile(x, 1)), max(x$mids))
}
