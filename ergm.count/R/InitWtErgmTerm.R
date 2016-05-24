#  File R/InitWtErgmTerm.R in package ergm.count, part of the Statnet suite
#  of packages for network analysis, http://statnet.org .
#
#  This software is distributed under the GPL-3 license.  It is free,
#  open source, and has the attribution requirements (GPL Section 7) at
#  http://statnet.org/attribution
#
#  Copyright 2008-2016 Statnet Commons
#######################################################################
InitWtErgmTerm.CMP<-function(nw, arglist, response, ...) {
  a <- check.ErgmTerm(nw, arglist,
                      varnames = NULL,
                      vartypes = NULL,
                      defaultvalues = list(),
                      required = NULL)
  list(name="sumlogfactorial",
       coef.names="CMP",
       inputs=NULL,
       dependence=FALSE,
       maxval=0)
}
