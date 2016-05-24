#  File R/control.stergm.R in package tergm, part of the Statnet suite
#  of packages for network analysis, http://statnet.org .
#
#  This software is distributed under the GPL-3 license.  It is free,
#  open source, and has the attribution requirements (GPL Section 7) at
#  http://statnet.org/attribution
#
#  Copyright 2003-2014 Statnet Commons
#######################################################################
###########################################################################
# The <control.stergm> function allows the ergm fitting process to be tuned
# by returning a list of several control parameters
############################################################################

control.stergm<-function(init.form=NULL,
                         init.diss=NULL,
                         init.method=NULL,
                         force.main = FALSE,                         

                         MCMC.prop.weights.form="default",MCMC.prop.args.form=NULL,
                         MCMC.prop.weights.diss="default",MCMC.prop.args.diss=NULL,
                         MCMC.init.maxedges=20000,
                         MCMC.init.maxchanges=20000,
                         MCMC.packagenames=c(),
                         
                         CMLE.MCMC.burnin = 1024*16,
                         CMLE.MCMC.interval = 1024,
                         CMLE.control=NULL,
                         CMLE.control.form=control.ergm(init=init.form, MCMC.burnin=CMLE.MCMC.burnin, MCMC.interval=CMLE.MCMC.interval, MCMC.prop.weights=MCMC.prop.weights.form, MCMC.prop.args=MCMC.prop.args.form, MCMC.init.maxedges=MCMC.init.maxedges, MCMC.packagenames=MCMC.packagenames, parallel=parallel, parallel.type=parallel.type, parallel.version.check=parallel.version.check, force.main=force.main),
                         CMLE.control.diss=control.ergm(init=init.diss, MCMC.burnin=CMLE.MCMC.burnin, MCMC.interval=CMLE.MCMC.interval, MCMC.prop.weights=MCMC.prop.weights.diss, MCMC.prop.args=MCMC.prop.args.diss, MCMC.init.maxedges=MCMC.init.maxedges, MCMC.packagenames=MCMC.packagenames, parallel=parallel, parallel.type=parallel.type, parallel.version.check=parallel.version.check, force.main=force.main),

                         CMLE.NA.impute=c(),
                         CMLE.term.check.override=FALSE,
                         
                         EGMME.main.method=c("Gradient-Descent"),

                         EGMME.MCMC.burnin.min=1000,
                         EGMME.MCMC.burnin.max=100000,
                         EGMME.MCMC.burnin.pval=0.5,
                         EGMME.MCMC.burnin.add=1,

                         MCMC.burnin=NULL, MCMC.burnin.mul=NULL,
                         
                         SAN.maxit=10,
                         SAN.control=control.san(coef=init.form,
                           SAN.prop.weights=MCMC.prop.weights.form,
                           SAN.prop.args=MCMC.prop.args.form,
                           SAN.init.maxedges=MCMC.init.maxedges,
                           
                           SAN.burnin=round(sqrt(EGMME.MCMC.burnin.min*EGMME.MCMC.burnin.max)),
                           SAN.packagenames=MCMC.packagenames,
                           
                           parallel=parallel,
                           parallel.type=parallel.type,
                           parallel.version.check=parallel.version.check),

                         SA.restarts=10,
                         
                         SA.burnin=1000,

                         # Plot the progress of the optimization.
                         SA.plot.progress=FALSE,
                         SA.max.plot.points=400,
                         SA.plot.stats=FALSE,
                         
                         # Initial gain --- if the process initially goes
                         # crazy beyond recovery, lower this.
                         SA.init.gain=0.1,
                         SA.gain.decay=0.5, # Gain decay factor.
                         
                         SA.runlength=25, # Number of jumps per .C call.

                         # Interval --- number of steps between
                         # successive jumps --- is computed
                         # adaptively.
                         SA.interval.mul=2, # Set the mean duration of extant ties this to be the interval.
                         SA.init.interval=500, # Starting interval.
                         SA.min.interval=20, # The lowest it can go.
                         SA.max.interval=500, # The highest it can go.


                         SA.phase1.minruns=4, # Number of pure-jitter runs before gradient calculation and such begin.
                         SA.phase1.tries=20, # Number of iterations of trying to find a reasonable configuration. FIXME: nothing happens if it's exceeded.
                         SA.phase1.jitter=0.1, # Initial jitter sd of each parameter..
                         SA.phase1.max.q=0.1, # FDR q-value for considering a gradient to be significant.
                         SA.phase1.backoff.rat=1.05, # If a run produces this relative increase in the objective function, it will be backed off.
                         SA.phase2.levels.max=40, # Maximum number of gain levels to go through.
                         SA.phase2.levels.min=4, # Minimum number of gain levels to go through.
                         SA.phase2.max.mc.se=0.001, # Maximum Monte-Carlo variation error of the parameter estimates as a fraction of total variation.
                         SA.phase2.repeats=400, # Maximum number of times gain a subphase can be repeated if the optimization is "going somewhere".
                         SA.stepdown.maxn=200, # Thin the draws for trend detection to get this many.
                         SA.stepdown.p=0.05, # If the combined p-value for the trend in the parameters is less than this, reset the subphase counter.
                         SA.stop.p=0.1, # If the combined p-value of the stopping tests exceeds this, finish the optimization.
                         SA.stepdown.ct=5, # Baseline number of times in a row the p-value must be above SA.stepdown.p to reduce gain.
                         SA.phase2.backoff.rat=1.1, # If a run produces this relative increase in the objective function, it will be backed off.
                         SA.keep.oh=0.5, # Fraction of optimization history that is used for gradient and covariance calculation.
                         SA.keep.min.runs=8, # Minimum number of runs that are used for gradient and covariance calculation.
                         SA.keep.min=0, # Minimum number of observations that are used for gradient and covariance calculation.
                         SA.phase2.jitter.mul=0.2, # The jitter standard deviation of each parameter is this times its standard deviation sans jitter.
                         SA.phase2.maxreljump=4, # Maximum jump per run, relative to the magnitude of other jumps in the history.
                         SA.guard.mul = 4, # The multiplier for the range of parameter values to compute the guard width.
                         SA.par.eff.pow = 1, # How do we scale rows of the estimating equation as a function of G?
                         SA.robust = FALSE, # Should the (slower) robust linear regression and covariance estimation be used?
                         SA.oh.memory = 100000, # Maximum number of jumps to store.
                         
                         SA.refine=c("mean","linear","none"), # Method, if any, used to refine the point estimate: linear interpolation, average, and none for the last value.
                         
                         SA.se=TRUE, # Whether to run Phase 3 to compute the standard errors.
                         SA.phase3.samplesize.runs=10, # This times the interval is the number of steps to estimate the standard errors.
                         SA.restart.on.err=TRUE, # Whether to wrap certain routines in try() statements so that they that an error is handled gracefully. Set to FALSE to debug errors in those routines.

                         seed=NULL,
                         parallel=0,
                         parallel.type=NULL,
                         parallel.version.check=TRUE){

  if(!is.null(MCMC.burnin) || !is.null(MCMC.burnin.mul)) stop("Control parameters MCMC.burnin and MCMC.burnin.mul are no longer used. See help for EGMME.MCMC.burnin.min, EGMME.MCMC.burnin.max, EGMME.MCMC.burnin.pval, EGMME.MCMC.burnin.pval, and CMLE.MCMC.burnin and CMLE.MCMC.interval for their replacements.")

  match.arg.pars=c("EGMME.main.method","SA.refine")

  if(!is.null(CMLE.control)) CMLE.control.form <- CMLE.control.diss <- CMLE.control
  
  control<-list()
  formal.args<-formals(sys.function())
  for(arg in names(formal.args))
    control[arg]<-list(get(arg))

  for(arg in match.arg.pars)
    control[arg]<-list(match.arg(control[[arg]][1],eval(formal.args[[arg]])))

  set.control.class()
}
