#  File R/control.simulate.stergm.R in package tergm, part of the Statnet suite
#  of packages for network analysis, http://statnet.org .
#
#  This software is distributed under the GPL-3 license.  It is free,
#  open source, and has the attribution requirements (GPL Section 7) at
#  http://statnet.org/attribution
#
#  Copyright 2003-2014 Statnet Commons
#######################################################################
########################################################################
# The <control.simulate.X> functions each create a list of paramaters
# for customizing simulation rountines
#
# --PARAMETERS--
#   prop.weights.form/
#   ptop.weighs.diss: specifies the method used to allocate probabilities
#                     of being proposed to dyads for the formation/dis-
#                     solution processes; options are "TNT",
#                     "random", and "default"; default="default" which
#                      picks a reasonable default considering any constraints
#   final           : whether only the final network of the simulation
#                     process should be returned (T or F); default= FALSE
#                     in which case, models, coefficients, stats matrices,
#                     and the toggle matrix are returned
#   maxchanges      : the maximum number of changes for which to allocate
#                     space; default=1000000
#
# --IGNORED--
#   prop.args.form/
#   prop.args.diss: an alternative, direct way of specifying additional
#                   arguments to proposal; as far as I can tell, the
#                   only use for 'prop.args' is to supply the name
#                   of a nodal attribute for use in the
#                   <InitMHP.nobetweengroupties> function, but this
#                   function will never be called in the path from
#                   <simulate.stergm> which is the only code using this
#                   control list.
#
# --RETURNED--
#   a list of the above parameters
#
#########################################################################

control.simulate.network<-function(MCMC.burnin.min=1000,
                                   MCMC.burnin.max=100000,
                                   MCMC.burnin.pval=0.5,
                                   MCMC.burnin.add=1,
                                   MCMC.burnin=NULL, MCMC.burnin.mul=NULL,
                                   MCMC.prop.weights.form="default",MCMC.prop.args.form=NULL,
                                   MCMC.prop.weights.diss="default",MCMC.prop.args.diss=NULL,                                  
                                   MCMC.init.maxedges=20000,
                                   MCMC.packagenames=c(),
                                   
                                   MCMC.init.maxchanges=1000000){
    if(!is.null(MCMC.burnin) || !is.null(MCMC.burnin.mul)) stop("Control parameters MCMC.burnin and MCMC.burnin.mul are no longer used. See help for EGMME.MCMC.burnin.min, EGMME.MCMC.burnin.max, EGMME.MCMC.burnin.pval, EGMME.MCMC.burnin.pval, and CMLE.MCMC.burnin and CMLE.MCMC.interval for their replacements.")
    control<-list()
    for(arg in names(formals(sys.function())))
      control[arg]<-list(get(arg))

    set.control.class()
  }

control.simulate.stergm<-function(MCMC.burnin.min=NULL,
                                  MCMC.burnin.max=NULL,
                                  MCMC.burnin.pval=NULL,
                                  MCMC.burnin.add=NULL,
                                  MCMC.burnin=NULL, MCMC.burnin.mul=NULL,
                                  MCMC.prop.weights.form=NULL,MCMC.prop.args.form=NULL,
                                  MCMC.prop.weights.diss=NULL,MCMC.prop.args.diss=NULL,                                  
                                  MCMC.init.maxedges=NULL,
                                  MCMC.packagenames=NULL,

                                  MCMC.init.maxchanges=NULL){
    if(!is.null(MCMC.burnin) || !is.null(MCMC.burnin.mul)) stop("Control parameters MCMC.burnin and MCMC.burnin.mul are no longer used. See help for EGMME.MCMC.burnin.min, EGMME.MCMC.burnin.max, EGMME.MCMC.burnin.pval, EGMME.MCMC.burnin.pval, and CMLE.MCMC.burnin and CMLE.MCMC.interval for their replacements.")
    control<-list()
    for(arg in names(formals(sys.function())))
      control[arg]<-list(get(arg))

    set.control.class()
  }


