#  File R/zzz.R in package ergm.userterms, part of the Statnet suite
#  of packages for network analysis, http://statnet.org .
#
#  This software is distributed under the GPL-3 license.  It is free,
#  open source, and has the attribution requirements (GPL Section 7) at
#  http://statnet.org/attribution
#
#  Copyright 2003-2013 Statnet Commons
#######################################################################
.onAttach <- function(lib, pkg){
  sm <- statnetStartupMessage("blkergm", c("statnet","ergm.userterm"), FALSE)
  if(!is.null(sm)){
    packageStartupMessage(sm)
    packageStartupMessage("An extension to the ergm package for fitting the block models based on the ergm.userterms")
  }
}



