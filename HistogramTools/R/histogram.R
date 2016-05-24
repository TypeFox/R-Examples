# Copyright 2011 Google Inc. All Rights Reserved.
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

.BreaksAreEquidistant <- function(breaks) {
  # Check if breaks are equally spaced.
  # If we had integer breakpoints we could just use
  # hist$equidist <- length(unique(diff(hist$breaks))) == 1
  diffs <- diff(breaks)
  all(abs(diffs - diffs[1]) < .Machine$double.eps^0.5 * max(diffs))
}

PreBinnedHistogram <- .BuildHistogram <- function(breaks, counts, xname="") {
  # Returns a histogram object from the given list of breaks and counts.

  stopifnot(is.numeric(breaks), is.numeric(counts))
  stopifnot(length(breaks) > 1)
  stopifnot(length(breaks) == (length(counts) + 1))
  hist <- list(breaks = breaks,
               counts = counts,
               density = counts / (sum(counts) * diff(breaks)),
               mids = (head(breaks, -1) + tail(breaks, -1)) / 2,
               xname = xname,
               equidist = .BreaksAreEquidistant(breaks))
  class(hist) <- "histogram"
  return(hist)
}

.UpdateHistogram <- function(x) {
  # Takes a histogram with possibly modified counts and breaks and updates
  # all other required fields to be consistent.

  stopifnot(inherits(x, "histogram"))
  stopifnot(is.numeric(x$breaks), is.numeric(x$counts))
  stopifnot(length(x$breaks) > 1)
  stopifnot(length(x$breaks) == (length(x$counts) + 1))
  x$density <- x$counts / (sum(x$counts) * diff(x$breaks))
  x$mids <-  (head(x$breaks, -1) + tail(x$breaks, -1)) / 2
  x$equidist <- .BreaksAreEquidistant(x$breaks)
  return(x)
}

.NewHistogramName <- function(x) {
  if(length(x) == 2) {
    return(paste("Merge of", x[[1]]$xname, "and", x[[2]]$xname))
  } else {
    paste("Merge of", length(x), "histograms")
  }
}

.AddManyHistograms <- function(x, main=.NewHistogramName(x)) {
  # Adds many histogram objects together that have the same bins.
  #
  # TODO(mstokely): Relax the requirement that the histograms have exactly
  # the same bins.
  #
  # Args:
  #   x: A list of S3 histogram objects
  #   main:  The name to set for the merged histogram (e.g. used in plots).
  #
  # Returns:
  #   An S3 histogram class suitable for plotting.
  stopifnot(all(sapply(x, inherits, "histogram")))
  br <- unname(lapply(x, function(y) y$breaks))
  if (!all((sapply(br, identical, y = br[[1]])))) {
    stop("Histograms must have identical breakpoints for this operation.")
  }
  # Now we know that all histograms have identical breaks.
  cnts <- unname(lapply(x, function(y) y$counts))
  sum.cnts <- Reduce("+", cnts)
  return(.BuildHistogram(breaks = x[[1]]$breaks,
                         counts = sum.cnts,
                         xname = main))
}

AddHistograms <- function(..., x=list(...), main=.NewHistogramName(x)) {
  # Adds two histogram objects that have the same bins.
  #
  # TODO(mstokely): Support histograms without exactly the same breaks.
  #
  # Args:
  #   ...: S3 histogram objects.
  #   x:   A list of S3 histogram objects.
  #   main:  The name to set for the merged histogram (e.g. used in plots).
  #
  # Returns:
  #   An S3 histogram class suitable for plotting.

  stopifnot(length(x) > 0)
  if (length(x) == 1) {
    return(x)
  } else {
    # length > 1
    return(.AddManyHistograms(x, main=main))
  }
}

ScaleHistogram <- function(x, factor=1 / Count(x)) {
  # Scale the counts of each bucket in a histogram by a factor.
  #
  # Args:
  #   x:       A histogram object.
  #   factor:  Scaling factor for the counts.
  #
  # Returns:
  #   An S3 histogram class suitable for plotting.

  stopifnot(inherits(x, "histogram"))
  stopifnot(is.numeric(factor), length(factor) == 1)
  x$counts <- x$counts * factor
  return(.UpdateHistogram(x))
}

# S3 Generics

as.histogram <- function(x, ...) {
  UseMethod("as.histogram")
}

# TODO(mstokely): RProtoBuf confusingly uses asMessage instead of as.Message
as.Message <- function(x) {
  UseMethod("as.Message")
}

as.histogram.Message <- function(x, ...) {
  # Converts a Protocol Buffer into an R Histogram.
  #
  # Args:
  #   x: An RProtoBuf Message of type HistogramTools.HistogramState.
  #
  # Returns:
  #   An S3 histogram class suitable for plotting.
  stopifnot(inherits(x, "Message"))

  if (x@type != "HistogramTools.HistogramState") {
    stop(paste("Unknown protocol message type", x@type, "only",
               "HistogramTools.HistogramState supported"))
  }
  hist <- .BuildHistogram(breaks = x$breaks, counts = x$counts)
  if (x$has("name")) {
    hist$xname <- x$name
  } else {
    hist$xname <- "HistogramTools.HistogramState"
  }
  return(hist)
}

as.Message.histogram <- function(x) {
  # Converts an R S3 histogram class into a HistogramState ProtoBuf.
  #
  # Args:
  #   x: An S3 histogram object.
  #
  # Returns:
  #   An RProtoBuf message of type HistogramTools.HistogramState
  stopifnot(inherits(x, "histogram"))
  # This catches NAs
  stopifnot(!is.null(x$breaks))
  stopifnot(is.numeric(x$breaks))
  if (!requireNamespace("RProtoBuf", quietly = TRUE)) {
      stop("RProtoBuf required to convert Histograms to Protocol Buffer Messages.")
  }

  # We can't conditionally require RProtoBuf and do this in onload()
  if (("RProtoBuf:DescriptorPool" %in% search()) &&
      !exists("HistogramTools.HistogramState",
              where="RProtoBuf:DescriptorPool")) {
    RProtoBuf::readProtoFiles(package="HistogramTools")
  }

  hist.class <- RProtoBuf::P("HistogramTools.HistogramState")
  # new is made generic in RProtoBuf, and we are not importing new in the
  # namespace so reference it explicitly here.
  hist.msg <- RProtoBuf::new(hist.class)

  hist.msg$counts <- x$counts
  hist.msg$breaks <- x$breaks
  hist.msg$name <- x$xname
  return(hist.msg)
}

# NB(mstokely): This causes an R CMD check warning about
# the histogram class being undocumented, so commented out.
#
# This would let us do :
#   as(hist(runif(100), plot=F), "Message")
#
# But we can already accomplish that with :
#   as.Message(hist(runif(100), plot=F))
#
# so it is not a big loss.
#setOldClass("histogram")
#setAs("histogram", "Message", as.Message.histogram)
