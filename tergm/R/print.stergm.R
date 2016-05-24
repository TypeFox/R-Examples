#  File R/print.stergm.R in package tergm, part of the Statnet suite
#  of packages for network analysis, http://statnet.org .
#
#  This software is distributed under the GPL-3 license.  It is free,
#  open source, and has the attribution requirements (GPL Section 7) at
#  http://statnet.org/attribution
#
#  Copyright 2003-2014 Statnet Commons
#######################################################################

print.stergm <- function (x, digits = max(3, getOption("digits") - 3), ...) {
  cat("Formation Coefficients:\n")
  print.default(format(x$formation.fit$coef, digits = digits), print.gap = 2, quote = FALSE)
  cat("Dissolution Coefficients:\n")
  print.default(format(x$dissolution.fit$coef, digits = digits), print.gap = 2, quote = FALSE)
}
