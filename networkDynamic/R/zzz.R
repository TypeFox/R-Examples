#  File networkDynamic/R/zzz.R
#  Part of the statnet package, http://statnetproject.org
#
#  This software is distributed under the GPL-3 license.  It is free,
#  open source, and has the attribution requirements (GPL Section 7) in
#    http://statnetproject.org/attribution
#
#  Copyright 2012 the statnet development team
######################################################################
# .onLoad is run when the package is loaded with library(networkDynamic)
#
######################################################################

.onAttach<- function(lib, pkg){
  sm <- statnetStartupMessage("networkDynamic",c('statnet'),TRUE)
  if(!is.null(sm)) packageStartupMessage(sm)
}
