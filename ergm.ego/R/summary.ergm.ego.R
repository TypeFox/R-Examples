#  File R/summary.ergm.ego.R in package ergm.ego, part of the Statnet suite
#  of packages for network analysis, http://statnet.org .
#
#  This software is distributed under the GPL-3 license.  It is free,
#  open source, and has the attribution requirements (GPL Section 7) at
#  http://statnet.org/attribution
#
#  Copyright 2015-2016 Statnet Commons
#######################################################################
## A thin wrapper around summary.ergm to get rid of a spurious error message.
.summary.ergm.ego <- function (object, ..., 
                          digits = max(3, getOption("digits") - 3),
                          correlation=FALSE, covariance=FALSE,
                          total.variation=TRUE){
  class(object) <- "ergm"
  summ <- NextMethod("summary")
  class(summ) <- c("summary.ergm.ego", "summary.ergm")
  summ
}

.print.summary.ergm.ego <- function (x, ...){
  print.summary.ergm(x, ..., print.deviances=FALSE)
}

