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

# TODO(mstokely): Consider: MergeAllEmptyBuckets() to compress empty buckets.

.MergeBucketsToBreakList <- function(x, breaks, FUN=sum) {
  # Merge adjacent buckets of a Histogram given a subset of the original breaks.
  #
  # This function combines adjacent buckets such that the returned
  # histogram has 'breaks' as breakpoints; breaks must be a subset of
  # the original breakpoints.
  #
  # Args:
  #   x:         An S3 histogram object
  #   breaks:    a vector giving the breakpoints between cells.
  #   FUN:       The function used to merge buckets.  Using anything other than
  #              sum here would break the density so use with care.
  #
  # Returns:
  #   An S3 histogram class suitable for plotting.

  stopifnot(is.numeric(breaks), length(breaks) > 1)
  if (!all(breaks %in% x$breaks)) {
    stop("List of breaks must be subset of existing breaks in histogram")
  }

  if (!isTRUE(all.equal(range(breaks), range(x$breaks)))) {
    stop("Range of new breakpoints too small.  Use SubsetHistogram first.")
  }
  # TODO(mstokely): Make this work and subset the histogram with a warning.
  #  if (max(breaks) < max(x$breaks)) {
  #    warning("Trimming buckets from histogram.")
  #    x <- SubsetHistogram(x, maxbreak=max(breaks))
  #  }
  #  if (min(breaks) < min(x$breaks)) {
  #    warning("Trimming buckets from histogram.")
  #    x <- SubsetHistogram(x, minbreak=min(breaks))
  #  }

  i <- which(x$breaks %in% breaks)
  bucket.grouping <- rep(head(breaks, -1), diff(i))
  tmp.df <- aggregate(x$counts, by=list(breaks=bucket.grouping), FUN)

  x$counts <- tmp.df$x
  x$breaks <- c(tmp.df$breaks, tail(x$breaks, 1))

  # Updated the other named list elements of the histogram class and return.
  return(.UpdateHistogram(x))
}

MergeBuckets <- function(x, adj.buckets=NULL, breaks=NULL, FUN=sum) {
  # Merge adjacent buckets of a Histogram.
  #
  # This only makes sense where the new bucket boundaries are a subset
  # of the previous bucket boundaries.  Only one of adj.buckets or
  # breaks should be specified.
  #
  # Args:
  #   x: An S3 histogram object
  #   adj.buckets: The number of adjacent buckets to merge.
  #   breaks: a vector giving the breakpoints between cells, or a
  #     single number giving number of cells.  Must have same range
  #     as x$breaks.
  #   FUN: The function used to merge buckets.  Using anything other than
  #     sum here would break the density so use with care.
  #
  # Returns:
  #   An S3 histogram class suitable for plotting.
  stopifnot(inherits(x, "histogram"))
  if (is.null(adj.buckets)) {
    stopifnot(is.numeric(breaks), length(breaks) > 0)
    if (length(breaks) > 1) {
      return(.MergeBucketsToBreakList(x, breaks, FUN))
    }
    stopifnot(breaks < length(x$breaks))
    # How many new buckets will we have.
    new.bucket.count <- breaks
    adj.buckets <- ceiling(length(x$counts) / new.bucket.count)
  } else {
    stopifnot(is.numeric(adj.buckets), length(adj.buckets) == 1,
              adj.buckets > 1)
    if (!is.null(breaks)) {
      stop("Only one of adj.buckets and breaks should be specified.")
    }
    new.bucket.count <- ceiling(length(x$counts) / adj.buckets)
  }

  # The last bucket may not be full, hence the length.out
  # TODO(mstokely): Consider bucket.grouping <- x$breaks[
  #                             ceiling(seq_along(x$counts) / adj.buckets)]
  # or: seq(from = 1, by = adj.buckets, length = new.bucket.count + 1)
  bucket.grouping <- rep(x$breaks[1+(0:new.bucket.count)*adj.buckets],
                         each=adj.buckets, length.out=length(x$counts))

  tmp.df <- aggregate(x$counts, by=list(breaks=bucket.grouping), FUN)

  x$counts <- tmp.df$x
  x$breaks <- c(tmp.df$breaks, tail(x$breaks, 1))

  # Updated the other named list elements of the histogram class and return.
  return(.UpdateHistogram(x))
}
