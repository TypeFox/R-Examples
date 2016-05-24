#  File R/stergm.CMLE.R in package tergm, part of the Statnet suite
#  of packages for network analysis, http://statnet.org .
#
#  This software is distributed under the GPL-3 license.  It is free,
#  open source, and has the attribution requirements (GPL Section 7) at
#  http://statnet.org/attribution
#
#  Copyright 2003-2014 Statnet Commons
#######################################################################
stergm.CMLE <- function(nw, formation, dissolution, constraints, times, offset.coef.form, offset.coef.diss,
                        eval.loglik,
                        estimate,
                        control,
                        verbose) {

  if(is.null(times)){
    # check if it is has network list class, or NOT a network and is a list (because all networks are lists)
    if(inherits(nw, "network.list") || (is.list(nw) & !is.network(nw))){
      times  <- seq_along(nw)
      warning("'times' argument was not provided to specify sampling time points for a list. Modeling transition between successive networks jointly. This behavior may change in the future.")
      if(var(sapply(nw,network.size))>0){
        stop("Networks in the network series passed must all be of the same size.")
      }
    }else if(inherits(nw,"networkDynamic")){
      times  <- c(0,1)
      warning("the 'times' argument was not provided to specify sampling time points for networkDynamic object. Modeling transition from time 0 to 1.")
      # TODO: should it check vertex activity to ensure that the effective sizes of the networks are the same at all time points
    }
  }
  
  if(length(times)<2) stop("Time points whose transition is to be modeled was not specified.")  


  # Construct a list of "from" networks and a list of "to" networks.
  if(inherits(nw,"networkDynamic")){

    # Grab only the needed vertices.
    subnw <- network.extract(nw, onset=min(times), terminus=max(times)+min(abs(diff(times)))/2)
    # Grab the vector of vertex activity indicators for each time
    # point, bind them into an n*T matrix. If any rows have
    # variability, we have a changing composition.
    if(any(apply(sapply(times, function(t) networkDynamic::is.active(subnw, at=t, v=1:network.size(subnw))), 1, var)>0)) warning("Active vertex set varies from time point to time point. Estimation may not work.")
    
    y0s <- lapply(times[-length(times)], function(t) networkDynamic::network.collapse(subnw, at=t, retain.all.vertices=TRUE))
    y1s <- lapply(times[-1], function(t) networkDynamic::network.collapse(subnw, at=t, retain.all.vertices=TRUE))
  }else if(inherits(nw, "network.list") || is.list(nw)){
    y0s <- nw[times[-length(times)]]
    y1s <- nw[times[-1]]
    
    if(!all(sapply(c(y0s,y1s),is.network))) stop("nw must be a networkDynamic, a network.list, or a list of networks.")
  }

  y0s <- lapply(y0s,standardize.network)
  y1s <- lapply(y1s,standardize.network)

  nwl <- c(y0s[1], y1s)
  class(nwl) <- "network.list"

  # Impute missing dyads in initial networks, if needed.
  y0s.NA <- sapply(y0s, network.naedgecount)>0 # Update which networks have missing dyads.
  y0s <- impute.network.list(y0s, control$CMLE.NA.impute, nwl.append=y1s[length(y1s)])
  y0s.NA <- sapply(y0s, network.naedgecount)>0 # Update which networks have missing dyads.
  
  if(any(y0s.NA)) stop("Transitioned-from network(s) at time(s) ", paste.and(times[-length(times)][y0s.NA]), " has missing dyads that cannot be imputed using the selected imputation options. Fix or add imputation options via CMLE.NA.impute control parameter.")

  if(length(times)>2){
    y0 <- combine.networks(y0s, standardized=TRUE,blockname=".stergm.CMLE.time.index")
    y1 <- combine.networks(y1s, standardized=TRUE,blockname=".stergm.CMLE.time.index")

    if(!control$CMLE.term.check.override){
      # Check that these networks can be combined for this model.

      eql <- function(target, current, tolerance = .Machine$double.eps^0.5)
          2*abs(target-current)/(abs(target)+abs(current))<tolerance # Some statistics don't give exact equality.
        
      bad.stat <-
        !(c(eql(apply(rbind(sapply(y0s, function(nw) summary(ergm.update.formula(formation,nw~., from.new="nw")))),1,sum), summary(ergm.update.formula(formation,y0~., from.new="y0"))),
            eql(apply(rbind(sapply(y0s, function(nw) summary(ergm.update.formula(dissolution,nw~., from.new="nw")))),1,sum), summary(ergm.update.formula(dissolution,y0~., from.new="y0")))) &
          c(eql(apply(rbind(sapply(y1s, function(nw) summary(ergm.update.formula(formation,nw~., from.new="nw")))),1,sum), summary(ergm.update.formula(formation,y1~., from.new="y1"))),
            eql(apply(rbind(sapply(y1s, function(nw) summary(ergm.update.formula(dissolution,nw~., from.new="nw")))),1,sum), summary(ergm.update.formula(dissolution,y1~., from.new="y1")))))
      if(any(bad.stat)) stop("Fitting the term(s) ", paste.and(unique(names(bad.stat)[bad.stat])), " over multiple network transitions is not supported at this time.")
    }
  }else{
    # We are about to do logical operations on networks, so make sure
    # tail-head orderings match up.
    y0 <- standardize.network(y0s[[1]])
    y1 <- standardize.network(y1s[[1]])
  }
  
  # Construct the formation and dissolution networks; the
  # network.update cannot be used to copy attributes from y0 to y.form and
  # y.diss, since NA edges will be lost.
  
  y.form <- y0 | y1
  y.form <- nvattr.copy.network(y.form, y0)
  formation <- ergm.update.formula(formation, y.form~., from.new="y.form")

  y.diss <- y0 & y1
  y.diss <- nvattr.copy.network(y.diss, y0)
  dissolution <- ergm.update.formula(dissolution, y.diss~., from.new="y.diss")

  # Construct new constraints

  constraints.form <- if(constraints==~.) ~atleast(y0) else ergm.update.formula(constraints, ~.+atleast(y0), from.new="y0")
  if(length(times)>2) constraints.form <- ergm.update.formula(constraints.form, ~.+blockdiag(".stergm.CMLE.time.index"))
  # TODO: Some unlucky variable names can break this. We need to figure out a way around this.
  environment(constraints.form) <- environment()
  
  constraints.diss <- if(constraints==~.) ~atmost(y0) else ergm.update.formula(constraints, ~.+atmost(y0), from.new="y0")
  if(length(times)>2) constraints.diss <- ergm.update.formula(constraints.diss, ~.+blockdiag(".stergm.CMLE.time.index"))
  # TODO: Some unlucky variable names can break this. We need to figure out a way around this.
  environment(constraints.diss) <- environment()
  
  # Apply initial values passed to control.stergm() the separate controls, if necessary.
  if(is.null(control$CMLE.control.form$init)) control$CMLE.control.form$init <- control$init.form
  if(is.null(control$CMLE.control.diss$init)) control$CMLE.control.diss$init <- control$init.diss
 
  model.form<-ergm.getmodel(formation, y.form, initialfit=TRUE)
  model.diss<-ergm.getmodel(dissolution, y.diss, initialfit=TRUE)

  if(!is.null(control$CMLE.control.form$init)){
    # Check length of control$CMLE.control.form$init.
    if(length(control$CMLE.control.form$init)!=length(model.form$etamap$offsettheta)) {
      if(verbose) cat("control$CMLE.control.form$init is", control$CMLE.control.form$init, "\n", "number of statistics is",length(model.form$coef.names), "\n")
      stop(paste("Invalid starting formation parameter vector control$CMLE.control.form$init:",
                 "wrong number of parameters."))
    }
  }else control$CMLE.control.form$init <- rep(NA, length(model.form$etamap$offsettheta)) # Set the default value of control$CMLE.control.form$init.
  if(!is.null(offset.coef.form)) control$CMLE.control.form$init[model.form$etamap$offsettheta]<-offset.coef.form
  names(control$CMLE.control.form$init) <- model.form$coef.names

  if(!is.null(control$CMLE.control.diss$init)){
    # Check length of control$CMLE.control.diss$init.
    if(length(control$CMLE.control.diss$init)!=length(model.diss$etamap$offsettheta)) {
      if(verbose) cat("control$CMLE.control.diss$init is", control$CMLE.control.diss$init, "\n", "number of statistics is",length(model.diss$coef.names), "\n")
      stop(paste("Invalid starting dissolution parameter vector control$CMLE.control.diss$init:",
                 "wrong number of parameters."))
    }
  }else control$CMLE.control.diss$init <- rep(NA, length(model.diss$etamap$offsettheta)) # Set the default value of control$CMLE.control.diss$init.  
  if(!is.null(offset.coef.diss)) control$CMLE.control.diss$init[model.diss$etamap$offsettheta]<-offset.coef.diss
  names(control$CMLE.control.diss$init) <- model.diss$coef.names


  # Translate the "estimate" from the stergm() argument to the ergm() argument.
  ergm.estimate <- switch(estimate,
                     CMLE = "MLE",
                     CMPLE = "MPLE")
  
  # Now, call the ergm()s:
  cat("Fitting formation...\n")
  fit.form <- ergm(formation, constraints=constraints.form, offset.coef=offset.coef.form, eval.loglik=eval.loglik, estimate=ergm.estimate, control=control$CMLE.control.form, verbose=verbose)
  cat("Fitting dissolution...\n")
  fit.diss <- ergm(dissolution, constraints=constraints.diss, offset.coef=offset.coef.diss, eval.loglik=eval.loglik, estimate=ergm.estimate, control=control$CMLE.control.diss, verbose=verbose)

  # Construct the output list. Conveniently, this is mainly a list consisting of two ergms.
  
  list(network=nwl, times=times, formation=formation, dissolution=dissolution, constraints=constraints, formation.fit=fit.form, dissolution.fit=fit.diss, estimate=estimate)
}
