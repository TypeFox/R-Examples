# Copyright 2012 Google Inc. All Rights Reserved.
# Author: mstokely@google.com (Murray Stokely)

TestTrimHistogram <- function() {
  hist.1 <- hist(c(1,2,3), breaks=0:9, plot=FALSE)
  hist.trimmed <- TrimHistogram(hist.1)
  # All non-zero elements still accounted for.
  checkTrue(sum(hist.1$counts) == sum(hist.trimmed$counts))
  # But fewer buckets:
  checkTrue(length(hist.1$counts) > length(hist.trimmed$counts))

  hist.2 <- hist(c(4,5,6), breaks=0:9, plot=FALSE)
  hist.trimmed <- TrimHistogram(hist.2)
  # All non-zero elements still accounted for.
  checkTrue(sum(hist.2$counts) == sum(hist.trimmed$counts))
  # But fewer buckets:
  checkTrue(length(hist.2$counts) > length(hist.trimmed$counts))

  # TODO(mstokely) This emits warnings.  Could check on that fact.
  zero.hist <- hist(numeric(), breaks=c(0,1,2,3,4,5,6,7,8,9), plot=FALSE)
  zero.trimmed <- TrimHistogram(zero.hist)
  # Don't do anything when all the buckets are empty.
  checkEquals(length(zero.hist$counts), length(zero.trimmed$counts))
}
