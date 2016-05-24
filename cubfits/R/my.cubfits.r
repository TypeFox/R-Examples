### Function to run MCMC inference for following model:
### Given phi, n, b; across n.G genes:
###   phi     ~ rlnorm(n.G, mu.Phi, sigma.Phi)
###   phi.Obs ~ rlnorm(n.G, log(phi), sigmaW)
###   y       ~ for(aa) rmultinom(n.G, invmlogit( phi * b[[aa]] ), n[[aa]] )
###
### Expects phi.Obs as vector of length n.G.
### 
### Subsequent Gibbs sampler:
### (1) Sample of b | phi, y, using VGAM fit for Gaussian proposal
### (2) Sample sigmaW | phi, phi.Obs
###     Sample (mu.Phi, sigma.Phi) | x 
### (3) Sample of phi | b, phi.Obs, nu.Phi, sigma.Phi, m.Phi, ww, sigmaW, y,
###     using logNormal proposal 

### All genes have observations.
my.cubfits <- function(reu13.df.obs, phi.Obs, y, n,
    nIter = 1000,
    b.Init = NULL, init.b.Scale = .CF.CONF$init.b.Scale,
        b.DrawScale = .CF.CONF$b.DrawScale,
        b.RInit = NULL,
    p.Init = NULL, p.nclass = .CF.CONF$p.nclass,
        p.DrawScale = .CF.CONF$p.DrawScale,
    phi.Init = NULL, init.phi.Scale = .CF.CONF$init.phi.Scale,
        phi.DrawScale = .CF.CONF$phi.DrawScale,
    model = .CF.CT$model[1],
    model.Phi = .CF.CT$model.Phi[1],
    adaptive = .CF.CT$adaptive[1],
    verbose = .CF.DP$verbose,
    iterThin = .CF.DP$iterThin,
    report = .CF.DP$report){

### Setup functions ###
  ### Setup function pointers by type or model.
  my.function <- my.init.function(model = model[1], model.Phi = model.Phi[1],
                                  adaptive = adaptive[1])
  my.ncoef <- my.function$my.ncoef

### Check Data ###
  ### check phi.Init
  if(is.null(phi.Init)){
    phi.Init <- phi.Obs + rlnorm(length(phi.Obs))
    phi.Init <- phi.Init / mean(phi.Init)
  }

  ### check phi.* is well-behaved
  my.check.data(phi.Obs = phi.Obs, phi.Init = phi.Init)

  ### Check if sort by ORF and length.
  my.check.rearrange(reu13.df.obs, y, n, phi.Obs = phi.Obs, phi.Init = phi.Init)

### Initial Storages ###
  ### Setup data structures for results
  n.G <- length(phi.Obs)                    # # of genes
  n.aa <- length(y)                         # # of amino acids
  nsyns <- sapply(y, function(ybit){ dim(ybit)[2] })
                                            # # of synomous codons
  nBparams <- my.ncoef * sum(nsyns - 1)     # total # of regression parameters
  nSave <- nIter / iterThin + 1             # # of space for iterations
  nPrior <- 3                               # # of prior parameters
  if(model.Phi == "logmixture"){
    nPrior <- 1 + 3 * p.nclass
  }
  if(.CF.CONF$estimate.bias.Phi){
    nPrior <- nPrior + 1                    # one more for bias.Phi
  }

  ### Storages for saving posterior samples.
  b.Mat <- my.generate.list(NA, nBparams, nSave)    # log(mu) and Delta.t
  p.Mat <- my.generate.list(NA, nPrior, nSave)      # prior parameters
  phi.Mat <- my.generate.list(NA, n.G, nSave)       # E[Phi | Phi^{obs}]
  logL.Mat <- my.generate.list(NA, 1, nSave)        # logL

### Initial Parameters ###
  ### Initial values for p first since scaling may change phi.Obs.
  p.Init <- my.pInit(p.Init, phi.Obs, model.Phi[1],
                     p.nclass = p.nclass, cub.method = "fits")

  ### Scaling.
  if(.CF.CONF$estimate.bias.Phi){
    phi.Obs.org <- phi.Obs
    phi.Obs <- phi.Obs / mean(phi.Obs)    # scale to mean 1
  }

  ### Initial values for b anb b.R.
  b.InitList <- .cubfitsEnv$my.fitMultinomAll(reu13.df.obs, phi.Init, y, n)
  if(is.null(b.RInit)){
    b.RInitList <- lapply(b.InitList, function(B){ B$R })
  } else{
    b.RInitList <- b.RInit
  }
  if(is.null(b.Init)){
    b.Init <- lapply(b.InitList,
                function(B){
                  B$coefficients +
                  init.b.Scale * backsolve(B$R, rnorm(nrow(B$R)))
                })
  } else{
    if(!is.null(b.Init[[1]]$R)){
      b.RInitList <- lapply(b.Init, function(B){ B$R })
    }
    b.Init <- lapply(b.Init, function(B){ B$coefficients })
  }
  b.InitVec <- unlist(b.Init)
  names(b.RInitList) <- names(reu13.df.obs)

  ### Initial values for training phi.
  if(is.null(phi.Init)){
    phi.Init <- phi.Obs * exp(init.phi.Scale * p.Init[1] * rnorm(n.G))
  }

  ### Push scaling back.
  if(.CF.CONF$estimate.bias.Phi){
    phi.Obs <- phi.Obs.org
  }

### Set current step ###
  ### Set current step for b.
  b.Mat[[1]] <- b.InitVec
  b.Curr <- b.Init

  ### Set current step for p.
  p.Mat[[1]] <- p.Init 
  p.Curr <- p.Init

  ### Set current step for phi.
  phi.Mat[[1]] <- phi.Init
  phi.Curr <- phi.Init
  log.phi.Obs <- log(phi.Obs)

  ### For hyper-prior parameters.
  hp.param <- list(log.phi.Obs.mean = mean(log.phi.Obs),
                   # hp.sigma.Phi = 1 / sqrt(var(log.phi.Obs)),
                   hp.Init = p.Init[-1])

  ### Set logL.
  logL.Curr <- -Inf
  if(.CF.CONF$compute.logL){
    tmp <- .cubfitsEnv$my.logLAll(phi.Curr, phi.Obs, y, n, b.Init,
                                  p.Curr, reu13.df = reu13.df.obs)
    logL.Curr <- sum(tmp)
    logL.Mat[[1]] <- logL.Curr 
  }

### MCMC here ###
  ### Get length for acceptance and adaptive storage.
  n.p <- 1
  if(.CF.CONF$estimate.bias.Phi){
    n.p <- n.p + 1
  }

  ### Set acceptance rate storage.
  my.set.acceptance(nIter + 1, n.aa, n.p = n.p, n.G = n.G)

  ### Set adaptive storage.
  if(.CF.CONF$estimate.bias.Phi && length(p.DrawScale) < n.p){
    ### Bias of phi is coupled with p parameters.
    p.DrawScale <- c(p.DrawScale, .CF.CONF$bias.Phi.DrawScale)
  }
  my.set.adaptive(nIter + 1,
                  n.aa = n.aa, b.DrawScale = b.DrawScale,
                  n.p = n.p, p.DrawScale = p.DrawScale,
                  n.G = n.G, phi.DrawScale = phi.DrawScale,
                  adaptive = adaptive[1])

  ### Run MCMC iterations
  my.verbose(verbose, 0, report)
  .cubfitsEnv$my.dump(0, list = c("b.Mat", "p.Mat", "phi.Mat", "logL.Mat"))

  ### MCMC start.
  for(iter in 1:nIter){
    ### Step 1: Update b using M-H step
    bUpdate <- .cubfitsEnv$my.drawBConditionalAll(
                 b.Curr, phi.Curr, y, n, reu13.df.obs,
                 b.RInitList = b.RInitList)
    b.Curr <- lapply(bUpdate, function(U){ U$bNew })

    ### Step 2: Draw other parameters.
    p.Curr <- .cubfitsEnv$my.pPropType(
                n.G, log.phi.Obs, phi.Curr, p.Curr, hp.param)

    ### Step 3: Update phi using M-H step.
    phi.Curr <- my.drawPhiConditionalAll(
                  phi.Curr, phi.Obs, y, n, b.Curr, p.Curr,
                  reu13.df = reu13.df.obs)

    ### Step logL:
    if(.CF.CONF$compute.logL && (iter %% iterThin) == 0){
      tmp <- .cubfitsEnv$my.logLAll(phi.Curr, phi.Obs, y, n, b.Curr,
                                    p.Curr, reu13.df = reu13.df.obs)
      logL.Curr <- sum(tmp)
    }

    ### Step A: Update scaling factor.
    if(iter %/% .CF.AC$renew.iter + 1 == .cubfitsEnv$curr.renew){
      my.copy.adaptive()
    } else{
      .cubfitsEnv$my.update.DrawScale(
        c("b", "p", "phi"),
        c(.CF.AC$b.DrawScale, .CF.AC$p.DrawScale, .CF.AC$phi.DrawScale))
    }

    ### Update Covariance Matrix
    if((iter/iterThin) %in% c(100, 200, 400, 800)){
      b.InitList <- .cubfitsEnv$my.fitMultinomAll(reu13.df.obs, phi.Curr, y, n)
      b.RInitList <- lapply(b.InitList, function(B){ B$R })
      names(b.RInitList) <- names(reu13.df.obs)
    }


    ### Dump parameters out.
    if((iter %% iterThin) == 0){
      thinnedIter <- iter / iterThin + 1
      b.Mat[[thinnedIter]] <- do.call("c", b.Curr)
      p.Mat[[thinnedIter]] <- p.Curr
      phi.Mat[[thinnedIter]] <- phi.Curr
      if(.CF.CONF$compute.logL){
        logL.Mat[[thinnedIter]] <- logL.Curr
      }
    }
    my.verbose(verbose, iter, report)
    .cubfitsEnv$my.dump(iter, list = c("b.Mat", "p.Mat", "phi.Mat", "logL.Mat"))
  } ### MCMC end.

### Check acceptance of last renew iteration.
  my.check.acceptance(c("b", "p", "phi"))

### Return ###
  aa.names <- names(y)
  in.names <- names(b.Mat[[1]])
  names(b.Mat[[1]]) <- mapBMatNames(in.names, aa.names, model = model)

  ret <- list(b.Mat = b.Mat, p.Mat = p.Mat, phi.Mat = phi.Mat,
              logL.Mat = logL.Mat,
              b.Init = b.Init, b.RInit = b.RInitList,
              p.Init = p.Init, phi.Init = phi.Init)
  ret
} # End of my.cubfits().

