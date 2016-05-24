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

ReadHistogramsFromDtraceOutputFile <- function(filename) {
  # Read the output of the dtrace tool and return a list of R histograms.
  #
  # The dtrace tool can output histograms of system quantities in a textual
  # representation.  This function parses those into R histogram objects for
  # further statistical analysis.
  #
  # Args:
  #   filename: The text output file of dtrace that includes distributions.
  #
  # Returns:
  #   A list of S3 histogram objects suitable for plotting.

  stopifnot(is.character(filename), length(filename) == 1)
  stopifnot(file.exists(filename))
  dtrace.text <- readLines(filename)
  dividers <- grepl("--- Distribution ---", dtrace.text)

  # lengths for each
  lengths <- diff(c(which(dividers), length(dtrace.text)))

  start <- which(dividers) - 1
  lengths1 <- lengths - 1
  myh <- sapply(1:length(start), function(i) {
    dtrace.text[start[i] + 0:lengths1[i]]
  })

  myhists <- lapply(myh, .BuildSingleHistogramFromDtraceOutput)
  names(myhists) <- lapply(myhists, function(x) str_trim(x$xname))
  return(myhists)
}

.BuildSingleHistogramFromDtraceOutput <- function(textlines) {
  # Build a single histogram from a portion of the text output of dtrace.
  #
  # The ReadHistogramsFromDtraceOutputFile() function breaks up the text
  # output of dtrace into individual chunks corresponding to different
  # distributions and then calls this function on each subset to generate
  # a single histogram.  textlines is of the form:
  #
  #[1] "  tcsh                                              "
  #[2] "           value  ------------- Distribution ------------- count    "
  #[3] "               0 |                                         0        "
  #[4] "               1 |@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ 4        "
  #[5] "               2 |                                         0        "
  #[6] ""
  #
  # Args:
  #   textlines: A character vector of the portion of dtrace output for one hist
  #
  # Returns:
  #   An S3 histogram object suitable for plotting.

  stopifnot(is.character(textlines), length(textlines) > 3)
  stopifnot(grepl("-- Distribution --", textlines[2]))
  textlines <- textlines[nzchar(textlines)]
  title <- textlines[1]
  headerline <- textlines[2]
  value.rightoffset <- regexpr("value", headerline, fixed=T)[1] + nchar("value")
  count.leftoffset <- regexpr("count", headerline, fixed=T)[1]
  bins <- as.numeric(sub("(^.*)\\|.*", "\\1", tail(textlines, -2)))
  counts <- as.numeric(sub("^.{59}(.*)", "\\1", tail(textlines, -2)))
  return(.BuildHistogram(breaks = bins,
                         # remove the last bin from counts, always 0.
                         counts = head(counts, -1),
                         xname = title))
}
