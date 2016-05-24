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

TestSubsetHistogram <- function() {
  hist.1 <- hist(c(1,1,2,2,7), breaks=0:9, plot=FALSE)
  hist.subset <- SubsetHistogram(hist.1, maxbreak=3)
  checkEquals(hist.subset$breaks, c(0, 1, 2, 3))
  checkEquals(hist.subset$counts, c(2, 2, 0))

  # Fail if we get get a break point out of range:
  checkException(SubsetHistogram(hist.1, minbreak=-100))
  checkException(SubsetHistogram(hist.1, maxbreak=100))

  # Or if its in range but not an existing breakpoint:
  checkException(SubsetHistogram(hist.1, minbreak=0.5))
  checkException(SubsetHistogram(hist.1, maxbreak=6.5))

  # Return original histogram if new breakpoints are existing ends.
  hist.nosubset <- SubsetHistogram(hist.1, minbreak=0)
  checkEquals(hist.nosubset$breaks, hist.1$breaks)
  checkEquals(hist.nosubset$counts, hist.1$counts)

  hist.nosubset <- SubsetHistogram(hist.1, maxbreak=9)
  checkEquals(hist.nosubset$breaks, hist.1$breaks)
  checkEquals(hist.nosubset$counts, hist.1$counts)
}
