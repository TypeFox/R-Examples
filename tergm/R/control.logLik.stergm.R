#  File R/control.logLik.stergm.R in package tergm, part of the Statnet suite
#  of packages for network analysis, http://statnet.org .
#
#  This software is distributed under the GPL-3 license.  It is free,
#  open source, and has the attribution requirements (GPL Section 7) at
#  http://statnet.org/attribution
#
#  Copyright 2003-2014 Statnet Commons
#######################################################################

control.logLik.stergm<-function(control.form = control.logLik.ergm(),
                                control.diss = control.logLik.ergm(),
                                control = NULL,
                                
                                seed=NULL){
  if(!is.null(control))
    control.form <- control.diss <- control
  
  control<-list()
  formal.args<-formals(sys.function())
  for(arg in names(formal.args))
    control[arg]<-list(get(arg))

  set.control.class()
}
