#  File R/InitReference.R in package ergm.count, part of the Statnet suite
#  of packages for network analysis, http://statnet.org .
#
#  This software is distributed under the GPL-3 license.  It is free,
#  open source, and has the attribution requirements (GPL Section 7) at
#  http://statnet.org/attribution
#
#  Copyright 2008-2016 Statnet Commons
#######################################################################
InitReference.Poisson <- function(lhs.nw, ...){
  list(name="Poisson")
}

InitReference.Binomial <- function(lhs.nw, trials, ...){
  
  list(name="Binomial", trials=trials)
}

InitReference.Geometric <- function(lhs.nw, ...){
  
  list(name="Geometric")
}
