#  File R/gof.ergm.ego.R in package ergm.ego, part of the Statnet suite
#  of packages for network analysis, http://statnet.org .
#
#  This software is distributed under the GPL-3 license.  It is free,
#  open source, and has the attribution requirements (GPL Section 7) at
#  http://statnet.org/attribution
#
#  Copyright 2015-2016 Statnet Commons
#######################################################################
# This file contains the following 8 functions for assessing goodness of fit
#         <gof>              <summary.gofobject>
#         <gof.default>      <plot.gofobject>
#         <gof.ergm>         <ergm.get.terms.formula>
#         <gof.formula>      <ergm.rhs.formula>
#=============================================================================



###############################################################################
# Each of the <gof.X> functions assesses the goodness of fit of X by comparison
# with 'control$nsim' ergm simulations of X
#
# --PARAMETERS--
#   object/formula: either an ergm object or a formula
#   ...           : additional parameters passed from within the program;
#                   these are ignored
#   init        : the parameters from which the simulations will be drawn;
#                   default=NULL;
#   control$nsim          : the number of simulated ergms, with which to compare X;
#                   default=100
#   burnin        : the number of proposals to disregard before any MCMC
#                   sampling is done; this is passed along to the simulation
#                   routines; default=10000
#   interval      : the number of proposals between sampled ergm statistics;
#                   this is passed along to the simulation rountines;
#                   default=1000
#   GOF           : a one-sided formula specifying which summary statistics
#                   should be used in the GOF comparison; choices include
#                       distance      espartners    dspartners
#                       odegree       idegree       degree
#                       triadcensus   model
#                   default=NULL; is internally mapped to 
#                   ~degree+espartners+distance if nw is undirected, and
#                   ~idegree+odegree+espartners+distance otherwise
#   constraints   : a one-sided formula of the constraint terms; options are
#                         bd        degrees        nodegrees
#                         edges     degreedist     idegreedist
#                         observed  odegreedist
#                   default="~ ."   
#   control       : a list of parameters for controlling GOF evaluation, as
#                   returned by <control.gof.X>; default=control.gof.X()
#                   (note that <control.gof.X> has different defaults 
#                    depending on the class of X)
#   seed          : an integer value at which to set the random generator;
#                   default=NULL
#   verbose       : whether to print information on the progress of the
#                   simulations; default=FALSE
#
# --RETURNED--
#   returnlist: a list with the following components for each term
#               G given in 'GOF'
#      summary.G: a matrix of summary statistics for the observed and
#                 simulated G's; if G takes on the values {G1, G2,...,Gq},
#                 the entries of 'summary.G' are
#         [i,1]-- the observed frequency of Gi
#         [i,2]-- the minimum value of Gi from the simulations
#         [i,3]-- the mean value of Gi from the simulations
#         [i,4]-- the maximum value of Gi from the simulations
#         [i,5]-- the p-value for the observed Gi estimated from the
#                 distribution of simulations
#      pobs.G   : a vector giving G's observed probability distribution
#      psim.G   : a matrix of G's simulated probability distributions; each
#                 row gives a distribution
#      bds.G    : the estimatd confidence interval, as the .025 and .975
#                 quantile values of the simulations
#      obs.G    : the vector of summary statistics for the observed X
#      sim.G    : the matrix of summary statistics for each simulated
#                 version of X
#
###############################################################################

gof.ergm.ego <- function (object, ..., 
                          GOF=c("model","degree"), 
                          control=control.gof.ergm(),
                          verbose=FALSE) {
  check.control.class(c("gof.ergm","gof.formula"))
  
  #Set up the defaults, if called with GOF==NULL
  if(is.null(GOF)){
    GOF<- ~degree
  }
  
  ## If a different constraint was specified, use it; otherwise, copy
  ## from the ERGM.

  control.transfer <- c("MCMC.burnin", "MCMC.interval", "MCMC.prop.weights", "MCMC.prop.args", "MCMC.packagenames", "MCMC.init.maxedges")
  for(arg in control.transfer)
    if(is.null(control[[arg]]))
      control[arg] <- list(object$control[[arg]])

  constraints <- object$constraints
  
  if(!is.null(control$seed)) {set.seed(as.integer(control$seed))}
  if (verbose) 
    cat("Starting GOF for the given ERGM formula.\n")

  GOF <- match.arg(GOF)
  
  # Calculate network statistics for the observed graph
  # Set up the output arrays of sim variables
  if(verbose)
    cat("Calculating observed network statistics.\n")

  n <- network.size(object$newnetwork)
  
  if ('model' %in% GOF) {
    obs.model <- summary(object$formula, scaleto=1)
    sim.model <- simulate(object, nsim=control$nsim, seed=control$seed, popsize=object$ppopsize, control=control.simulate.ergm.ego(simulate.control=set.control.class("control.simulate.formula",control)),...,verbose=verbose, statsonly=TRUE)/n
  }

  if ('degree' %in% GOF) {
    egodata <- object$egodata
    s <- 
    obs.deg <- summary(as.formula(paste0("egodata~degree(0:",n-1,")")), scaleto=1)
    sim.deg <- simulate(object, nsim=control$nsim, seed=control$seed, popsize=object$ppopsize, control=control.simulate.ergm.ego(simulate.control=set.control.class("control.simulate.formula",control)),...,verbose=verbose, statsonly=TRUE, monitor=as.formula(paste0("~degree(0:",n-1,")")))
    sim.deg <- sim.deg[, ncol(sim.deg)-((n-1):0), drop=FALSE]/n
  }

  if(verbose)
    cat("Starting simulations.\n")

  # calculate p-values
  
  returnlist <- list(network.size=n, GOF=as.formula(paste0("~",GOF)))
  
  if ('model' %in% GOF) {
    pval.model <- apply(sim.model <= obs.model[col(sim.model)],2,mean)
    pval.model.top <- apply(sim.model >= obs.model[col(sim.model)],2,mean)
    pval.model <- cbind(obs.model,apply(sim.model, 2,min), apply(sim.model, 2,mean),
                        apply(sim.model, 2,max), pmin(1,2*pmin(pval.model,pval.model.top)))
    dimnames(pval.model)[[2]] <- c("obs","min","mean","max","MC p-value")
    pobs.model <- pval.model.top
    psim.model <- apply(sim.model,2,rank)/nrow(sim.model)
    bds.model <- apply(psim.model,2,quantile,probs=c(0.025,0.975))

    returnlist$summary.model <- returnlist$pval.model <- pval.model
    returnlist$pobs.model <- pobs.model
    returnlist$psim.model <- psim.model
    returnlist$bds.model <- bds.model
    returnlist$obs.model <- obs.model
    returnlist$sim.model <- sim.model
  }

  if ('degree' %in% GOF) {
    pval.deg <- apply(sim.deg <= obs.deg[col(sim.deg)],2,mean)
    pval.deg.top <- apply(sim.deg >= obs.deg[col(sim.deg)],2,mean)
    pval.deg <- cbind(obs.deg,apply(sim.deg, 2,min), apply(sim.deg, 2,mean),
                      apply(sim.deg, 2,max), pmin(1,2*pmin(pval.deg,pval.deg.top)))
    dimnames(pval.deg)[[2]] <- c("obs","min","mean","max","MC p-value")
    pobs.deg <- obs.deg/sum(obs.deg)
    psim.deg <- sweep(sim.deg,1,apply(sim.deg,1,sum),"/")
    psim.deg[is.na(psim.deg)] <- 1
    bds.deg <- apply(psim.deg,2,quantile,probs=c(0.025,0.975))
    
    returnlist$summary.deg <- returnlist$pval.deg <- pval.deg
    returnlist$pobs.deg <- pobs.deg
    returnlist$psim.deg <- psim.deg
    returnlist$bds.deg <- bds.deg
    returnlist$obs.deg <- obs.deg
    returnlist$sim.deg <- sim.deg
  }

  class(returnlist) <- "gofobject"
  returnlist
}
