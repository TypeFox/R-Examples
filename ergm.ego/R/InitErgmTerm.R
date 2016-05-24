#  File R/InitErgmTerm.R in package ergm.ego, part of the Statnet suite
#  of packages for network analysis, http://statnet.org .
#
#  This software is distributed under the GPL-3 license.  It is free,
#  open source, and has the attribution requirements (GPL Section 7) at
#  http://statnet.org/attribution
#
#  Copyright 2015-2016 Statnet Commons
#######################################################################
# This is just an alias for the number of edges, to distinguish it in the summary tables.

InitErgmTerm.netsize.adj<-function(nw, arglist, ...) {
  a <- check.ErgmTerm(nw, arglist,
                      varnames = NULL,
                      vartypes = NULL,
                      defaultvalues = list(),
                      required = NULL)
  
  list(name="edges", coef.names="netsize.adj", dependence=FALSE,
       minval = 0, maxval = network.dyadcount(nw,FALSE), conflicts.constraints="edges", pkgname="ergm")
}
