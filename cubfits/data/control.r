### cubfits control variables.

### For method controls.
.CF.CT <- list(
  model = c("roc"),                             # main models
  type.p = c("lognormal_RW",
             "lognormal_fix",
             "lognormal_bias",
             "logmixture"),                     # proposal for hyperparameters
  type.Phi = c("RW_Norm"),                      # proposal for Phi
  model.Phi = c("lognormal", "logmixture"),     # prior of Phi
  init.Phi = c("PM"),                           # initial methods for Phi
  init.fit = c("RW_Norm", "current", "random"), # how is b proposed
  parallel = c("lapply", "mclapply",
               "task.pull", "pbdLapply"),       # parallel functions
  adaptive = c("simple", "none"),                # method for adaptive mcmc
  prior.dist.M = c("uniform", "normal"),
  prior.dist.Sphi = c("uniform", "normal")
)

### For configuration of initial and draw scaling.
.CF.CONF <- list(
  scale.phi.Obs = FALSE,           # if phi.Obs were scaled to mean 1
  init.b.Scale = 1,                # initial b scale
  init.phi.Scale = 1,              # initial phi scale
  p.nclass = 2,                    # # of classes if mixture phi
  b.DrawScale = 1,                 # drawing scale for b if random walk
  p.DrawScale = 0.1,               # drawing scale for p if random walk
  phi.DrawScale = 1,               # random walk scale for phi
  phi.pred.DrawScale = 1,          # random walk scale for phi.pred
  sigma.Phi.DrawScale = 1,         # random walk scale for sigma.Phi
  bias.Phi.DrawScale = 0.1,        # random walk scale for bias.Phi
  estimate.bias.Phi = FALSE,        # if estimate bias of phi during MCMC
  estimate.Phi.noise = TRUE,       # estimate the noise in the phi observed data (sigma epsilon)
  estimate.S.Phi = TRUE,           # if FALSE, sPhi is fixed accross the run
  compute.logL = TRUE              # if compute logL in each iteration
)

### For optimization.
.CF.OP <- list(
  optim.method = c("Brent"),                       # for optim()
  stable.min.exp = .Machine$double.max.exp * 0.1,  # minimum exponent
  stable.max.exp = .Machine$double.max.exp * 0.5,  # maximum exponent
  lower.optim = 1e-4,                              # lower of d logL(x)
  upper.optim = 1e2,                               # upper of d logL(x)
  lower.integrate = 0.0,                           # lower of \int L(x)
  upper.integrate = Inf                            # upper of \int L(x)
)

### For dumpping data.
.CF.DP <- list(
  dump = FALSE,                    # if dumping within MCMC
  iter = 1000,                     # iterations per dumping
  prefix.dump = "dump_",           # path and file names of dumping
  verbose = FALSE,                 # if verbose
  iterThin = 1,                    # iterations to thin chain
  report = 10,                     # iterations to report
  report.proc = 100                # iterations to report proc.time()
)

### For addaptive control.
.CF.AC <- list(
  renew.iter = 100,                # per renewing iterations
  target.accept.lower = 0.2,       # target acceptant rate lower bound
  target.accept.upper = 0.3,       # target acceptant rate upper bound
  scale.increase = 1.2,            # increase scale size
  scale.decrease = 0.8,            # decrease scale size
  sigma.lower = 1e-2,              # lower bound of relative scale size
  sigma.upper = 1e2                # upper bound of relative scale size
)

### For parameters as reestimated for Yeast according to Yassour's data.
.CF.PARAM <- list(
  # phi.meanlog = -0.441473,         # yassour mean for log(phi)
  # phi.sdlog = 1.393285,            # yassour sd for log(phi)
  phi.meanlog = -1.125,            # mean of log(phi), -s^2/2
  phi.sdlog = 1.5,                  # sd of log(phi)
  prior.M.a = 0,                    # first parameter of density function of prior (e.g. dnorm(x, mean=a, sd=b))
  prior.M.b = 0.35,                     # second parameter of density function of prior (e.g dnorm(x, mean=a, sd=b))
  prior.Sphi.a = 0.52,          # meanlog for lognormal Sphi prior
  prior.Sphi.b = 0.33          # sdlog for lognormal Sphi prior
)
