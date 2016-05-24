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

.TrimRightZeroBuckets <- function(x) {
  # Trims the empty buckets in the right-tail of a histogram.
  #
  # Args:
  #   x: An S3 histogram object.
  # Returns:
  #   A new S3 histogram object with any empty buckets in the
  #   right-tail removed.

  # Copy over max/min/sum of squares/class, etc.
  # Override breaks, counts, and midpoints to trim off the rightmost
  # zero buckets.

  if (all(x$counts == 0)) {
    warning("All buckets of histogram zero, returning unmodified.")
    return(x)
  }
  biggest.nonzero <- max(which(x$counts > 0))
  if (biggest.nonzero == length(x$counts)) {
    # Last bucket is non-zero.
    return(x)
  }
  # Since histogram objects are named lists, we can use within
  # to override the necessary elements.
  new.hist <- within(unclass(x), {
    breaks <- head(breaks, biggest.nonzero + 1)
    counts <- head(counts, biggest.nonzero)
    mids <- head(mids, biggest.nonzero)
    density <- head(density, biggest.nonzero)
  })
  attributes(new.hist) <- attributes(x)
  return(new.hist)
}

.TrimLeftZeroBuckets <- function(x) {
  # Trims the empty buckets in the left-tail of a histogram.
  #
  # Args:
  #   x: An S3 histogram object.
  # Returns:
  #   A new S3 histogram object with any empty buckets in the
  #   left-tail removed.

  if (all(x$counts == 0)) {
    warning("All buckets of histogram zero, returning unmodified.")
    return(x)
  }
  smallest.nonzero <- which.max(x$counts > 0)
  if (smallest.nonzero == 1) {
    # First bucket is non-zero.
    return(x)
  }
  new.hist <- within(unclass(x), {
    breaks <- tail(breaks, length(breaks) - smallest.nonzero + 1)
    counts <- tail(counts, length(counts) - smallest.nonzero + 1)
    mids <- tail(mids, length(mids) - smallest.nonzero + 1)
    density <- tail(density, length(density) - smallest.nonzero + 1)
  })
  attributes(new.hist) <- attributes(x)
  return(new.hist)
}

TrimHistogram <- function(x, left=TRUE, right=TRUE) {
  # Trims empty buckets from the left and right tails of a histogram.
  #
  # Args:
  #   x:     An R histogram object.
  #   left:  boolean, if TRUE, remove empty buckets from left tail of histogram.
  #   right: boolean, if TRUE, remove empty buckets from right tail of histogram.
  # Returns:
  #   x:  An R histogram object, possibly with fewer buckets.
  stopifnot(inherits(x, "histogram"))
  tmp.hist <- x
  if (left) {
    tmp.hist <- .TrimLeftZeroBuckets(tmp.hist)
  }
  if (right) {
    tmp.hist <- .TrimRightZeroBuckets(tmp.hist)
  }
  return(tmp.hist)
}
