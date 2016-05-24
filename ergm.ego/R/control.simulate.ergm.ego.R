#  File R/control.simulate.ergm.ego.R in package ergm.ego, part of the Statnet suite
#  of packages for network analysis, http://statnet.org .
#
#  This software is distributed under the GPL-3 license.  It is free,
#  open source, and has the attribution requirements (GPL Section 7) at
#  http://statnet.org/attribution
#
#  Copyright 2015-2016 Statnet Commons
#######################################################################
control.simulate.ergm.ego <- function(
  ppop.wt = c("round","sample"),
  SAN.control = control.san(),
  simulate.control = control.simulate(),
  ...){
  match.arg.pars <- c("ppop.wt")

  control<-list()
  formal.args<-formals(sys.function())
  formal.args[["..."]]<-NULL
  for(arg in names(formal.args))
    control[arg]<-list(get(arg))
  if(length(list(...))) stop("Unrecognized control parameter: ",arg,".")

  for(arg in match.arg.pars)
    control[arg]<-list(match.arg(control[[arg]][1],eval(formal.args[[arg]])))

  set.control.class("control.simulate.ergm.ego")
}
