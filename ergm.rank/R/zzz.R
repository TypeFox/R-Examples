#  File R/zzz.R in package ergm.rank, part of the Statnet suite
#  of packages for network analysis, http://statnet.org .
#
#  This software is distributed under the GPL-3 license.  It is free,
#  open source, and has the attribution requirements (GPL Section 7) at
#  http://statnet.org/attribution
#
#  Copyright 2008-2016 Statnet Commons
#######################################################################
.onAttach <- function(lib, pkg){
  sm <- statnetStartupMessage("ergm.rank",c("statnet"),FALSE)
  if(!is.null(sm)) packageStartupMessage(sm)
}

.onLoad <- function(lib, pkg){
  .RegisterMHPs()
  .RegisterConstraintImplications()
  .RegisterInitMethods()
}

.RegisterMHPs <- function(){
  ergm.MHP.table("c", "CompleteOrder", "",  0, "random", "AlterSwap")
}

.RegisterConstraintImplications <- function(){
}

.RegisterInitMethods <- function(){
  ergm.init.methods("CompleteOrder", c("CD","zeros"))
}
