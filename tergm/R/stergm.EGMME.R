#  File R/stergm.EGMME.R in package tergm, part of the Statnet suite
#  of packages for network analysis, http://statnet.org .
#
#  This software is distributed under the GPL-3 license.  It is free,
#  open source, and has the attribution requirements (GPL Section 7) at
#  http://statnet.org/attribution
#
#  Copyright 2003-2014 Statnet Commons
#######################################################################
################################################################################
# The <stergm> function fits stergms from a specified formation and dissolution
# formula returning approximate MLE's based on MCMC estimation.
#
# --PARAMETERS--
#   formation   : the formation formula, as 'nw ~ term(s)'
#   dissolution : the dissolution formula, as 'nw ~ term(s)'
#   theta.form0 : the intial theta formation coefficients, or optionally if
#                 these are to be estimates, the string "MPLE";
#                 default="MPLE"
#   theta.diss  : the initial theta dissolution coefficients
#   seed        : an integer starting value for the random number generator;
#                 default=NULL
#   MH.burnin   : the number of proposals used in each MCMC step; this is ignored
#                 unless 'control$main.method'="Robbins-Monro"; any other style or
#                 the default style will not recognize this parameter;
#                 default=1000
#   constraints : a one-sided formula of the constraint terms; options are
#                      bd        degrees        nodegrees
#                      edges     degreedist     indegreedist
#                      observed  outdegreedist
#                default="~ ."; these may not work currently.
#   target.stats   :  a vector of the mean value parameters;
#                  default=the observed statistic from the 'nw' in formula
#   control     :  a list of control parameters returned from <control.stergm>;
#                  default=<control.stergm>()
#   verbose     :  whether ergm should be verbose (T or F); default=FALSE
#
# --RETURNED--
#   because a stergm object is the return type of several functions, and
#   because this is a rather lengthy list, and because the returned items
#   of this function borrow from the other stergm.* functions, this list 
#   provides the returned items for all funtions returning a stergm.
#   The symbol preceding each component indicates which function returns it,
#   but remember that, <stergm> will additionally return the items from
#   one of the other stergm functions as well:
#       <stergm>             = $
#       <stergm.RM>          = @
#
#   the components include:
#
#     @   coef.form   : the estimated formation model coefficients
#     @   coef.diss   : the estimated dissolution model coefficients
#      &  eta         : the estimated formation ?? coefficients
#   $     offset      : a logical vector whose ith entry tells whether the
#                        ith curved theta coeffient was offset/fixed
#   $     etamap      :  the list constituting the theta->eta mapping for the
#                        formation model; for details of its components,
#                        see <ergm.etamap>
#   $     MH.burnin   :  the number of proposals made in each MCMC step
#   $     formation   : the formation formula, as 'nw ~ term(s)'
#   $     dissolution : the dissolution formula, as 'nw ~ term(s)'
#   $     constraints : the constraints formula
#     @&  newnetwork  :  the final network sampled
#     @&  network    :  the 'nw' inputted to <ergm> via the 'formula'
#   $     prop.args.form     :  the MHP formation arguments passed to the
#                               InitMHP rountines
#   $     prop.args.diss     :  the MHP dissolution arguments passed to the
#                               InitMHP rountines
#   $     prop.weights.form  :  the method used to allocate probabilities of
#                               being proposed to dyads in the formation stage,
#                               as "TNT", "random", "nonobserved", or "default"
#      &  theta.original     :  the theta values at the start of the MCMC 
#                               sampling
#     @   theta.form.original:  the formation theta values at the start of the
#                               MCMC sampling
#   $     prop.weights.diss  :  as 'prop.weights.form', but for the dissolution
#                               model
#
################################################################################

stergm.EGMME <- function(nw, formation, dissolution, constraints, offset.coef.form, offset.coef.diss,
                   targets, target.stats, estimate,
                 control,
                 verbose) {

  if(!is.network(nw)) stop("Argument nw must be a network.")

#  if(is.null(nw %n% "lasttoggle")) nw %n% "lasttoggle" <- rep(round(-.Machine$integer.max/2), network.dyadcount(nw))
	
  if(is.null(nw %n% "time")) nw %n% "time" <- 0

  # EGMME requires targets, or there will be an error
  if (is.null(targets)) stop('stergm.EGMME requires targets parameter be specified')
 
  # Allow the user to specify targets as copied from formation or dissolution formula.
  if(is.character(targets)){
    targets <- switch(targets,
                      formation = formation,
                      dissolution = dissolution)
  }
  
  if(length(targets)==3){
    warning("Targets formula has an LHS, which will be ignored in favor of nw.")
    targets <- targets[c(1,3)] # in a formula f<-y~x, f[1]=~, f[2]=y, and f[3]=x
  }

  targets <- ergm.update.formula(targets,nw~., from.new="nw")
  formation <- ergm.update.formula(formation,nw~., from.new="nw")
  dissolution <- ergm.update.formula(dissolution,nw~., from.new="nw")

  # target formula should not have offsets. removing them
  if (any(offset.info.formula(targets)$term)) {
    message("Targets formula should not contain offset terms;
                they have been been removed.")
    targets <- remove.offset.formula(targets)
  }

  control.transfer <- list(EGMME.MCMC.burnin.min="MCMC.burnin.min",
                           EGMME.MCMC.burnin.max="MCMC.burnin.max",
                           EGMME.MCMC.burnin.pval="MCMC.burnin.pval",
                           EGMME.MCMC.burnin.add="MCMC.burnin.add")
  for(arg in names(control.transfer))
      if(is.null(control[[control.transfer[[arg]]]]))
          control[control.transfer[[arg]]] <- list(control[[arg]])

  
  if(!is.null(target.stats)){

    nw.stats<-summary(targets)
    if(length(nw.stats)!=length(target.stats))
      stop("Incorrect length of the target.stats vector: should be ", length(nw.stats), " but is ",length(target.stats),".")
    
    if (!is.null(control$init.form) && length(target.stats) != length(control$init.form)) {
      if (length(target.stats) != length(control$SAN.control$coef)) {
        cat('SAN initial coefficients should have same length as targets. Setting them to SAN defaults.\n')
        control$SAN.control$coef <- NULL
      }
    }
        
    if(verbose) cat("Constructing an approximate response network.\n")
    ## If target.stats are given, overwrite the given network and targets
    ## with SAN-ed network and targets.
    
    newnw <- try({
      for(srun in seq_len(control$SAN.maxit)){
        nw<-suppressWarnings(
          san(targets, target.stats=target.stats,
                control=control$SAN.control,
                constraints=constraints,
                verbose=verbose))
        targets<-ergm.update.formula(targets,nw~., from.new="nw")
        nw.stats <- summary(targets)
        srun <- srun + 1
        if(verbose){
          cat(paste("Finished SAN run",srun,"\n"))
        }
        if(verbose){
          cat("SAN summary statistics:\n")
          print(nw.stats)
          cat("Meanstats Goal:\n")
          print(target.stats)
          cat("Difference: SAN target.stats - Goal target.stats =\n")
          print(round(nw.stats-target.stats,0))
        }
        if(sum((nw.stats-target.stats)^2) <= 5) break
      }
      nw
    }, silent=TRUE)
    if(inherits(newnw,"try-error")){
      cat("SAN failed or is not applicable. Increase burn-in if there are problems.\n")
    }else nw <- newnw
    formation <- ergm.update.formula(formation,nw~., from.new="nw")
    dissolution <- ergm.update.formula(dissolution,nw~., from.new="nw")
  }

  if (verbose) cat("Initializing Metropolis-Hastings proposals.\n")
  MHproposal.form <- MHproposal(constraints, weights=control$MCMC.prop.weights.form, control$MCMC.prop.args.form, nw, class="f")
  MHproposal.diss <- MHproposal(constraints, weights=control$MCMC.prop.weights.diss, control$MCMC.prop.args.diss, nw, class="d")

  if(!is.dyad.independent(MHproposal.form$arguments$constraints) || !is.dyad.independent(MHproposal.diss$arguments$constraints)){
    warning("Dyad-dependent constraint imposed. Note that the constraint is applied to the post-formation and post-dissolution networks y+ and y-, not the next time-step's network. This behavior may change in the future.")
  }
  
  model.form <- ergm.getmodel(formation, nw, expanded=TRUE, role="formation")
  model.diss <- ergm.getmodel(dissolution, nw, expanded=TRUE, role="dissolution")
  model.mon <- ergm.getmodel(targets, nw, expanded=TRUE, role="target")

  if(any(model.form$etamap$canonical==0) || any(model.diss$etamap$canonical==0) || any(model.mon$etamap$canonical==0)) stop("Equilibrium GMME for models based on curved ERGMs is not supported at this time.")

  p.free<-sum(!model.form$etamap$offsettheta)+sum(!model.diss$etamap$offsettheta)
  if(p.free==0) stop("Model specification has no free parameters (all are offsets).")
  q<-length(model.mon$etamap$offsettheta)
  if(p.free>q) stop("Fitting ",p.free," free parameters on ",q," target statistics. The specification is underidentified.")

  model.mon$nw.stats <- summary(model.mon$formula)
  model.mon$target.stats <- if(!is.null(target.stats)) vector.namesmatch(target.stats, names(model.mon$nw.stats)) else model.mon$nw.stats

  # If some control$init is specified...
  
  if(!is.null(control$init.form)){
    # Check length of control$init.form.
    if(length(control$init.form)!=length(model.form$etamap$offsettheta)) {
      if(verbose) cat("control$init.form is", control$init.form, "\n", "number of statistics is",length(model.form$coef.names), "\n")
      stop(paste("Invalid starting formation parameter vector control$init.form:",
                 "wrong number of parameters."))
    }
  }else control$init.form <- rep(NA, length(model.form$etamap$offsettheta)) # Set the default value of control$init.form.
  if(!is.null(offset.coef.form)) control$init.form[model.form$etamap$offsettheta]<-offset.coef.form
  names(control$init.form) <- model.form$coef.names

  if(!is.null(control$init.diss)){
    # Check length of control$init.diss.
    if(length(control$init.diss)!=length(model.diss$etamap$offsettheta)) {
      if(verbose) cat("control$init.diss is", control$init.diss, "\n", "number of statistics is",length(model.diss$coef.names), "\n")
      stop(paste("Invalid starting dissolution parameter vector control$init.diss:",
                 "wrong number of parameters."))
    }
  }else control$init.diss <- rep(NA, length(model.diss$etamap$offsettheta)) # Set the default value of control$init.diss.  
  if(!is.null(offset.coef.diss)) control$init.diss[model.diss$etamap$offsettheta]<-offset.coef.diss
  names(control$init.diss) <- model.diss$coef.names

  initialfit <- stergm.EGMME.initialfit(control$init.form, control$init.diss, nw, model.form, model.diss, model.mon, control, verbose)
  
  if(verbose) cat("Fitting STERGM Equilibrium GMME.\n")

  if(control$parallel){
    cl <- ergm.getCluster(control, verbose=verbose)
    if(verbose && !is.null(cl)) cat("Using parallel cluster.\n")
    on.exit(suppressWarnings(try(ergm.stopCluster(cl),silent=TRUE)))
  }else cl <- NULL
  
  Cout <- switch(control$EGMME.main.method,
                 "Gradient-Descent" = stergm.EGMME.GD(initialfit$formation.fit$coef,
                   initialfit$dissolution.fit$coef, nw, model.form, model.diss, model.mon,
                   control=control, MHproposal.form=MHproposal.form,
                   MHproposal.diss=MHproposal.diss, cl=cl,
                  verbose),
                 stop("Method ", control$EGMME.main.method, " is not implemented.")
                )

  if(!is.null(cl)) ergm.stopCluster(cl)

  out <- list(network = nw, formation = formation, dissolution = dissolution, targets = targets, target.stats=model.mon$target.stats, estimate=estimate, covar = Cout$covar, opt.history=Cout$opt.history, sample=Cout$sample, sample.obs=NULL, control=control, reference = ~Bernoulli, mc.se = Cout$mc.se, constraints = constraints,
              formation.fit = with(Cout, list(network=nw, formula=formation, coef = eta.form, covar=covar.form, etamap = model.form$etamap, offset = model.form$etamap$offsettheta, constraints=constraints, estimate=estimate, control=control, reference = ~Bernoulli, mc.se = mc.se.form)),
              dissolution.fit = with(Cout, list(network=nw, formula=dissolution, coef = eta.diss, covar=covar.diss, etamap = model.diss$etamap, offset = model.diss$etamap$offsettheta, constraints=constraints, estimate=estimate, control=control, reference = ~Bernoulli, mc.se = mc.se.diss))
              )
  class(out$formation.fit)<-class(out$dissolution.fit)<-"ergm"
  
  out
}
