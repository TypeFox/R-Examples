# Copyright 2014 Google Inc. All rights reserved.
#
# Use of this source code is governed by a BSD-style
# license that can be found in the LICENSE file or at
# http://opensource.org/licenses/BSD-3-Clause

IfElse <- function(test, yes, no, ...) {
  # Conditional evaluation, like 'test ? yes : no' in C, or
  # 'if(test) yes else no' in R.
  #
  # Args:
  #   test: a logical value. Missing values are not supported.
  #   yes: any object; this is returned if test is TRUE
  #   no:  any object; this is returned if test is FALSE;
  #        or, a logical value, if there are ... arguments.
  #   ... additional arguments; there should be an even number of them, e.g.
  #       IfElse(test1, yes1, test2, yes2, no2)
  #
  # For comparison,
  #   ifelse(test = logical vector, yes = a vector, no = a vector)
  #     The result is a vector with values from the second and third arguments.
  #   IfElse(logical scalar, yes = object, no = object)
  #     The result is either the second or third argument, in its entirety.
  #     Only one of yes and no is evaluated.
  if (test) {
    yes
  } else if (missing(..1)) {
    no
  } else {
    Recall(test = no, ...)
  }
}

if (FALSE) { # manual testing code
  IfElse(TRUE, 1:2, 3:4)
  IfElse(FALSE, 1:2, 3:4)
  IfElse(NA, 1:2, 3:4) # not supported, fails
  IfElse(FALSE, 1:2, TRUE, 3:5, 6:8)
  IfElse(FALSE, 1:2, FALSE, 3:5, 6:8)
  IfElse(FALSE, 1:2, FALSE, 3:5, TRUE, 6:9, 10:13)
  IfElse(FALSE, 1:2, FALSE, 3:5, FALSE, 6:9, 10:13)
}
