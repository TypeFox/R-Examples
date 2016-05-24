#  File R/stergm.EGMME.initialfit.R in package tergm, part of the Statnet suite
#  of packages for network analysis, http://statnet.org .
#
#  This software is distributed under the GPL-3 license.  It is free,
#  open source, and has the attribution requirements (GPL Section 7) at
#  http://statnet.org/attribution
#
#  Copyright 2003-2014 Statnet Commons
#######################################################################
stergm.EGMME.initialfit<-function(init.form, init.diss, nw, model.form, model.diss, model.mon, control, verbose=FALSE){
  if(!is.null(control$init.method) && control$init.method == "zeros"){
    init.form[is.na(init.form)]<-0
    init.diss[is.na(init.diss)]<-0
  }else if(!any(is.na(init.form)) && !any(is.na(init.diss))){
    # Don't need to do anything.
  }else if(all(model.form$coef.names[!model.form$etamap$offsettheta] %in% model.mon$coef.names)
           && (
                all(model.diss$etamap$offsettheta)
                || (
                     model.diss$coef.names == "edges"
                     && "mean.age" %in% model.mon$coef.names
                    )
                )
           && all(model.diss$coef.names %in% model.form$coef.names)
           && is.dyad.independent(model.diss$formula)){
    if(verbose) cat("Formation statistics are analogous to targeted statistics, dissolution is fixed or is edges with a mean.age target, dissolution terms appear to have formation analogs, and dissolution process is dyad-independent, so using edges dissolution approximation  (Carnegie et al.).\n")

    if(!all(model.diss$etamap$offsettheta)){ # This must mean that the two provisos above are satisfied.
      mean.age <- model.mon$target.stats[model.mon$coef.names=="mean.age"]
      init.diss <- log(mean.age+1)
      names(init.diss) <- "edges"
    }
    
    # Fit an ERGM to the formation terms:
    form.targets <- model.mon$target.stats[match(model.form$coef.names,model.mon$coef.names)]
    form.targets <- form.targets[!model.form$etamap$offsettheta]
    init.form<-coef(ergm(model.form$formula,control=control.ergm(init=init.form), target.stats=form.targets, eval.loglik=FALSE))
    # Now, match up non-offset formation terms with dissolution terms.
    # In case it's not obvious (it's not to me) what the following
    # does, it takes non-offset elements of init.form, then, from
    # those, it takes those elements that correspond (by name) to the
    # dissolution coefficient names and decrements them by init.diss.
    #
    # Yes, I am also a little surprised that assigning to a
    # double-index works.
    init.form[!model.form$etamap$offsettheta][match(names(init.diss),names(init.form[!model.form$etamap$offsettheta]))] <-
      init.form[!model.form$etamap$offsettheta][match(names(init.diss),names(init.form[!model.form$etamap$offsettheta]))] - init.diss
  }else{
    stop("No initial parameter method for specified model and targets combination is implemented. Specify via control$init.form and control$init.diss .")
  }
  out <- list(formation = model.form$formula, dissolution = model.diss$formula, targets = model.mon$formula, target.stats=model.mon$target.stats, nw = nw, control = control, formation.fit = list(coef=init.form, etamap = model.form$etamap), dissolution.fit = list(coef=init.diss, etamap = model.diss$etamap))
  class(out)<-"stergm"
  out
}
