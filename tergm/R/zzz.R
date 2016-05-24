#  File R/zzz.R in package tergm, part of the Statnet suite
#  of packages for network analysis, http://statnet.org .
#
#  This software is distributed under the GPL-3 license.  It is free,
#  open source, and has the attribution requirements (GPL Section 7) at
#  http://statnet.org/attribution
#
#  Copyright 2003-2014 Statnet Commons
#######################################################################
.onAttach <- function(lib, pkg){
  sm <- statnetStartupMessage("tergm", c("statnet"), FALSE)
  if(!is.null(sm)) packageStartupMessage(sm)
}

.onLoad <- function(lib, pkg){
  .RegisterMHPs()
  .RegisterConstraintImplications()
}

.RegisterMHPs <- function(){
  ergm.MHP.table("c", "Bernoulli", "atleast",  0, "random", "formationMLE")
  ergm.MHP.table("c", "Bernoulli", "atleast+bd",  0, "random", "formationMLE")
  ergm.MHP.table("c", "Bernoulli", "atleast",  1, "TNT", "formationMLETNT")
  ergm.MHP.table("c", "Bernoulli", "atleast+bd",  1, "TNT", "formationMLETNT")

  ergm.MHP.table("c", "Bernoulli", "atmost",  0, "random", "dissolutionMLE")
  ergm.MHP.table("c", "Bernoulli", "atmost+bd",  0, "random", "dissolutionMLE")
  ergm.MHP.table("c", "Bernoulli", "atmost",  1, "TNT", "dissolutionMLETNT")
  ergm.MHP.table("c", "Bernoulli", "atmost+bd",  1, "TNT", "dissolutionMLETNT")

  ergm.MHP.table("c", "Bernoulli", "atleast+observed",  0, "random", "formationNonObservedMLE")
  ergm.MHP.table("c", "Bernoulli", "atleast+bd+observed",  0, "random", "formationNonObservedMLE")
  ergm.MHP.table("c", "Bernoulli", "atleast+observed",  1, "TNT", "formationNonObservedMLETNT")
  ergm.MHP.table("c", "Bernoulli", "atleast+bd+observed",  1, "TNT", "formationNonObservedMLETNT")
  ergm.MHP.table("c", "Bernoulli", "atmost+observed",  0, "random", "dissolutionNonObservedMLE")
  ergm.MHP.table("c", "Bernoulli", "atmost+bd+observed",  0, "random", "dissolutionNonObservedMLE")
  ergm.MHP.table("c", "Bernoulli", "atmost+observed",  1, "TNT", "dissolutionNonObservedMLETNT")
  ergm.MHP.table("c", "Bernoulli", "atmost+bd+observed",  1, "TNT", "dissolutionNonObservedMLETNT")

  ergm.MHP.table("c", "Bernoulli", "atleast+blockdiag",  0, "random", "formationMLEblockdiag")
  ergm.MHP.table("c", "Bernoulli", "atleast+bd+blockdiag",  0, "random", "formationMLEblockdiag")
  ergm.MHP.table("c", "Bernoulli", "atleast+blockdiag",  1, "TNT", "formationMLEblockdiagTNT")
  ergm.MHP.table("c", "Bernoulli", "atleast+bd+blockdiag",  1, "TNT", "formationMLEblockdiagTNT")
  ergm.MHP.table("c", "Bernoulli", "atmost+blockdiag",  0, "random", "dissolutionMLEblockdiag")
  ergm.MHP.table("c", "Bernoulli", "atmost+bd+blockdiag",  0, "random", "dissolutionMLEblockdiag")
  ergm.MHP.table("c", "Bernoulli", "atmost+blockdiag",  1, "TNT", "dissolutionMLEblockdiagTNT")
  ergm.MHP.table("c", "Bernoulli", "atmost+bd+blockdiag",  1, "TNT", "dissolutionMLEblockdiagTNT")

  ergm.MHP.table("c", "Bernoulli", "atleast+blockdiag+observed",  0, "random", "formationNonObservedMLEblockdiag")
  ergm.MHP.table("c", "Bernoulli", "atleast+bd+blockdiag+observed",  0, "random", "formationNonObservedMLEblockdiag")
  ergm.MHP.table("c", "Bernoulli", "atleast+blockdiag+observed",  1, "TNT", "formationNonObservedMLEblockdiagTNT")
  ergm.MHP.table("c", "Bernoulli", "atleast+bd+blockdiag+observed",  1, "TNT", "formationNonObservedMLEblockdiagTNT")
  ergm.MHP.table("c", "Bernoulli", "atmost+blockdiag+observed",  0, "random", "dissolutionNonObservedMLEblockdiag")
  ergm.MHP.table("c", "Bernoulli", "atmost+bd+blockdiag+observed",  0, "random", "dissolutionNonObservedMLEblockdiag")
  ergm.MHP.table("c", "Bernoulli", "atmost+blockdiag+observed",  1, "TNT", "dissolutionNonObservedMLEblockdiagTNT")
  ergm.MHP.table("c", "Bernoulli", "atmost+bd+blockdiag+observed",  1, "TNT", "dissolutionNonObservedMLEblockdiagTNT")
  
  ergm.MHP.table("f", "Bernoulli", "",  0, "random", "formation")
  ergm.MHP.table("f", "Bernoulli", "bd",  0, "random", "formation")
  ergm.MHP.table("f", "Bernoulli", "",  1, "TNT", "formationTNT")
  ergm.MHP.table("f", "Bernoulli", "bd",  1, "TNT", "formationTNT")
  ergm.MHP.table("d", "Bernoulli", "",  0, "random", "dissolution")
  ergm.MHP.table("d", "Bernoulli", "bd",  0, "random", "dissolution")
  ergm.MHP.table("d", "Bernoulli", "",  1, "TNT", "dissolutionTNT")
  ergm.MHP.table("d", "Bernoulli", "bd",  1, "TNT", "dissolutionTNT")

}

.RegisterConstraintImplications <- function(){
  ergm.ConstraintImplications("atleast", c())
  ergm.ConstraintImplications("atmost", c())
}
