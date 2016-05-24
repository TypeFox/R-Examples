#  File R/stergm.R in package tergm, part of the Statnet suite
#  of packages for network analysis, http://statnet.org .
#
#  This software is distributed under the GPL-3 license.  It is free,
#  open source, and has the attribution requirements (GPL Section 7) at
#  http://statnet.org/attribution
#
#  Copyright 2003-2014 Statnet Commons
#######################################################################
################################################################################
# stergm --- fit Separable Temporal ERGMs.
################################################################################

stergm <- function(nw, formation, dissolution, constraints = ~., estimate, times=NULL, offset.coef.form=NULL, offset.coef.diss=NULL,
                   targets=NULL, target.stats=NULL,
                   eval.loglik=TRUE,
                   control=control.stergm(),
                   verbose=FALSE, ...) {
  check.control.class()
  
  if(!is.null(control$seed))  set.seed(as.integer(control$seed))

  estimate <- match.arg(estimate,c("CMLE","CMPLE","EGMME"))
  
  if(!inherits(formation,"formula") || !inherits(dissolution,"formula"))
    stop("Arguments formation and dissolution must be formulas.")

  if(length(formation)==3){
    warning("Formation formula has an LHS, which will be ignored in favor of nw.")
    formation <- formation[c(1,3)] # in a formula f<-y~x, f[1]=~, f[2]=y, and f[3]=x
  }

  if(length(dissolution)==3){
    warning("Dissolution formula has an LHS, which will be ignored in favor of nw.")
    dissolution <- dissolution[c(1,3)] # in a formula f<-y~x, f[1]=~, f[2]=y, and f[3]=x
  }
  
  # lasttoggle
  if(estimate=="EGMME"){
  duration.dependent <- is.lasttoggle(nw,formation,dissolution,targets)
  
  if(duration.dependent)
    nw %n% "lasttoggle" <- NVL(nw %n% "lasttoggle",rep(round(-.Machine$integer.max/2), network.dyadcount(nw)))  else nw %n% "lasttoggle" <- NULL
  }
  
  out <- switch(estimate,
                CMLE=,
                CMPLE=stergm.CMLE(nw, formation, dissolution, constraints,
                  times, offset.coef.form, offset.coef.diss, eval.loglik,
                  estimate, control, verbose),
                EGMME=stergm.EGMME(nw, formation, dissolution, constraints,
                  offset.coef.form, offset.coef.diss,
                  targets, target.stats, estimate, control, verbose)
                  )
  
  
  out$formation <- formation
  out$dissolution <- dissolution
  out$control <- control
  
  class(out)<-"stergm"
  out
}
