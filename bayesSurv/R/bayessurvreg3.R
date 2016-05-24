################################################
#### AUTHOR:     Arnost Komarek             ####
####             (2005)                     ####
####                                        ####
#### FILE:       bayessurvreg3.R            ####
####                                        ####
#### FUNCTIONS:  bayessurvreg3              ####
################################################

### ======================================
### bayessurvreg3
### ======================================
## Survival regression with G-spline random intercept
## and G-spline error
##
## 01/02/2005: start working on it
## 02/02/2005: finished
## 07/12/2006: extension allowing inclulsion of the estimated correlation between the onset and time-to-event random intercepts
## 31/05/2013: extension allowing for misclassification of the event status added
##             (cooperation with Maria Jose Garcia Zattera and Alejandro Jara)
##
bayessurvreg3 <- function
(  formula,
   random,
   formula2,
   random2,
   data = parent.frame(),
   classification,
   classParam = list(Model = c("Examiner", "Factor:Examiner"),
                     a.sens = 1, b.sens = 1, a.spec = 1, b.spec = 1,
                     init.sens = NULL, init.spec = NULL),
   na.action = na.fail,
   onlyX = FALSE,
   nsimul = list(niter = 10, nthin = 1, nburn = 0, nwrite = 10),   
   prior,
   prior.beta,
   prior.b,
   init = list(iter = 0),
   mcmc.par = list(type.update.a = "slice",   k.overrelax.a = 1,   k.overrelax.sigma = 1,   k.overrelax.scale = 1,
                   type.update.a.b = "slice", k.overrelax.a.b = 1, k.overrelax.sigma.b = 1, k.overrelax.scale.b = 1),
   prior2,
   prior.beta2,
   prior.b2,
   init2,
   mcmc.par2 = list(type.update.a = "slice",   k.overrelax.a = 1,   k.overrelax.sigma = 1,   k.overrelax.scale = 1,
                    type.update.a.b = "slice", k.overrelax.a.b = 1, k.overrelax.sigma.b = 1, k.overrelax.scale.b = 1),
   priorinit.Nb,
   rho = list(type.update = "fixed.zero", init = 0, sigmaL = 0.1),
   store = list(a = FALSE, a2 = FALSE, y = FALSE, y2 = FALSE, r = FALSE, r2 = FALSE, b = FALSE, b2 = FALSE,
                a.b = FALSE, a.b2 = FALSE, r.b = FALSE, r.b2 = FALSE), 
   dir = getwd())
{  
  thispackage = "bayesSurv"
  #thispackage = NULL
  store  <- bayessurvreg3.checkStore(store)
  nsimul <- bayessurvreg.checknsimul(nsimul)
  
  transform = function(t){log(t)}
  dtransform = function(t){1/t}
  transform2 = function(t){t}
  dtransform2 = function(t){1}
  
  ## Give a function call to be recorded in a resulting object.
  call <- match.call(expand.dots = TRUE)
  
  ## Extract all the design information from the function call
  ## Do not transform the response at this moment
  doubly <- ifelse(missing(formula2), FALSE, TRUE)
  m <- match.call(expand.dots = FALSE)
  des <- bayessurvreg.design(m = m, formula = formula, random = random, data = data, transform = transform2, dtransform = transform2)
  if (doubly) des2 <- bayessurvreg.design(m = m, formula = formula2, random = random2, data = data, transform = transform2, dtransform = transform2)
  else        des2 <- list(nX = 0, n = des$n, nrandom = 0, randomInt = FALSE)
  
  if (onlyX){
    if (doubly) return(list(X = des$X, X2 = des2$X))
    else        return(des$X)
  } 

  if (!des$nrandom){
    store$b   <- FALSE
    store$a.b <- FALSE
    store$r.b <- FALSE
  }
  if (!des2$nrandom){
    store$b2   <- FALSE
    store$a.b2 <- FALSE
    store$r.b2 <- FALSE
  }  

  ## Perform some checks
  if (des$nrandom > 1) stop("Only intercept may appear in 'random'")
  if (des$nrandom & !des$randomInt) stop("Only intercept may appear in 'random'")
  
  if (des2$nrandom > 1) stop("Only intercept may appear in 'random2'")
  if (des2$nrandom & !des2$randomInt) stop("Only intercept may appear in 'random2'")

  nobs <- des$n
  if (doubly){
    nobs2 <- des2$n
    if (nobs != nobs2) stop("Inconsistent formula and formula2 (different number of observations indicated)")
  }

  
  ### ------------------------------------------------------------------------------------
  ### 201305: Code related to misclassification models
  ### ------------------------------------------------------------------------------------
  if (missing(classification)){
    mclass <- list(Model = "NONE", nModel = 0, nvisit = rep(1, nobs),
                   nExaminer = 0, labelExaminer = "0", Examiner = rep(0, nobs),
                   nFactor   = 0, labelFactor   = "0", Factor   = rep(0, nobs),
                   Prior     = rep(1, 4),
                   logTime   = 0, Status = 0,
                   init.sens = 0, init.spec = 0)

    ## mclass$nModel == 0 ==> no misclassification model 
    
  }else{
    if (doubly) stop("misclassification model not (yet) implemented for doubly-interval-censored data.")    
    mclass <- list()        ## To keep the information being then passed to C++
    
    if (!is.data.frame(classification)) stop("classification must be a data.frame")
    if (ncol(classification) < 4) stop("too low number of columns in the classification data.frame")    
    if (ncol(classification) == 4){
      colnames(classification) <- c("idUnit", "Time", "Examiner", "Status")
      classification[, "Factor"] <- 0
    }else{
      classification <- classification[, 1:5]
      colnames(classification) <- c("idUnit", "Time", "Examiner", "Status", "Factor")
    }  

    ### Some checks
    if (any(is.na(classification[, "idUnit"])))          stop("NA's not allowed in the idUnit column of classification.")    
    if (any(is.na(classification[, "Time"])))            stop("NA's not allowed in the Time column of classification.")
    if (any(classification[, "Time"] <= 0))              stop("all values in the Time column of classification must be positive.")
    if (any(is.na(classification[, "Examiner"])))        stop("NA's not allowed in the Examiner column of classification.")
    if (any(is.na(classification[, "Status"])))          stop("NA's not allowed in the Status column of classification.")
    if (any(!(classification[, "Status"] %in% c(0, 1)))) stop("all values in the Status column of classification must be 0/1.")
    
    ### idUnit values that appear in the data
    IDUnit <- unique(classification[, "idUnit"])                                             ### This maintains original ordering

    ### Consistency check    
    if (length(IDUnit) != nobs) stop("classification data.frame not consistent with the survival data (different number of experimental units)")

    ### Misclassification model
    classParam$Model <- classParam$Model[1]
    mclass$Model <- classParam$Model
    mclass$nModel <- pmatch(mclass$Model, table = c("Examiner", "Factor:Examiner"))
    if (is.na(mclass$nModel)) stop("unknown classParam$Model specified.")
        
    ### Number of visits per unit (tooth)
    mclass$nvisit <- as.numeric(table(factor(classification[, "idUnit"], levels = IDUnit)))  ### Also here, the original ordering is maintained

    ### Number of examinators
    mclass$nExaminer <- length(unique(classification[, "Examiner"]))

    ### Examiners labeled 0, ..., nExaminer - 1 and their original labels
    mclass$labelExaminer <- levels(factor(classification[, "Examiner"]))
    mclass$Examiner <- as.numeric(factor(classification[, "Examiner"])) - 1

    ### Factor related variables (if needed)
    switch (mclass$Model,
      "Examiner" = {
        mclass$nFactor <- 1

        mclass$labelFactor <- "0"
        mclass$Factor <- rep(0, sum(mclass$nvisit))
      },
      "Factor:Examiner" = {
        ### Number of unique values of Factor
        mclass$nFactor <- length(unique(classification[, "Factor"]))

        ### Factor values labeled 0, ..., nFactor - 1 and their original labels
        mclass$labelFactor <- levels(factor(classification[, "Factor"]))
        mclass$Factor <- as.numeric(factor(classification[, "Factor"])) - 1            
      },
      stop("some part of the code not yet implemented for this option.")             
    )
    
    ### Prior parameters
    if (is.null(classParam$a.sens)) classParam$a.sens <- 1
    if (is.null(classParam$b.sens)) classParam$b.sens <- 1    
    if (is.null(classParam$a.spec)) classParam$a.spec <- 1
    if (is.null(classParam$b.spec)) classParam$b.spec <- 1
    if (classParam$a.sens <= 0 | classParam$b.sens <= 0 | classParam$a.spec <= 0 | classParam$b.spec <= 0) stop("incorrect prior parameter for misclassification sensitivities/specificities")
    mclass$Prior <- c(classParam$a.sens, classParam$b.sens, classParam$a.spec, classParam$b.spec)
    names(mclass$Prior) <- c("a.sens", "b.sens", "a.spec", "b.spec")    

    ### Transformation of the visit times
    mclass$logTime <- log(classification[, "Time"])

    ### Status
    mclass$Status <- classification[, "Status"]
    
    ### Initial sensitivities
    if (is.null(classParam$init.sens)){
      switch (mclass$Model,
        "Examiner" = {
          classParam$init.sens <- runif(mclass$nExaminer, min = 0.80, max = 0.90)
          names(classParam$init.sens) <- mclass$labelExaminer
        },
        "Factor:Examiner" = {
          classParam$init.sens <- matrix(runif(mclass$nFactor * mclass$nExaminer, min = 0.80, max = 0.90), nrow = mclass$nFactor, ncol = mclass$nExaminer)
          rownames(classParam$init.sens) <- mclass$labelFactor
          colnames(classParam$init.sens) <- mclass$labelExaminer
        },
        stop("some part of the code not yet implemented for this option.") 
      )
    }
    if (length(classParam$init.sens) != mclass$nFactor * mclass$nExaminer) stop("incorrect classParam$init.sens specified")
    if (any(is.na(classParam$init.sens))) stop("NA's not allowed in classParam$init.sens")
    if (any(classParam$init.sens <= 0) | any(classParam$init.sens >= 1)) stop("all classParam$init.sens must lie in (0, 1)")

    mclass$init.sens <- as.numeric(classParam$init.sens)
    
    ### Initial specificities
    if (is.null(classParam$init.spec)){
      switch (mclass$Model,
        "Examiner" = {
          classParam$init.spec <- runif(mclass$nExaminer, min = 0.80, max = 0.90)
          names(classParam$init.spec) <- mclass$labelExaminer
        },
        "Factor:Examiner" = {
          classParam$init.spec <- matrix(runif(mclass$nFactor * mclass$nExaminer, min = 0.80, max = 0.90), nrow = mclass$nFactor, ncol = mclass$nExaminer)
          rownames(classParam$init.spec) <- mclass$labelFactor
          colnames(classParam$init.spec) <- mclass$labelExaminer
        },
        stop("some part of the code not yet implemented for this option.") 
      )
    }        
    if (length(classParam$init.spec) != mclass$nFactor * mclass$nExaminer) stop("incorrect classParam$init.spec specified")
    if (any(is.na(classParam$init.spec))) stop("NA's not allowed in classParam$init.spec")
    if (any(classParam$init.spec <= 0) | any(classParam$init.spec >= 1)) stop("all classParam$init.spec must lie in (0, 1)")

    mclass$init.spec <- as.numeric(classParam$init.spec)    
  }  
  ### ------------------------------------------------------------------------------------
    
  ## Priors and inits for beta parameters
  if (missing(init))        init <- list()
  if (missing(init2))       init2 <- list()
  if (missing(mcmc.par))    mcmc.par <- list()
  if (missing(mcmc.par2))   mcmc.par2 <- list()
  if (!doubly)              mcmc.par2 <- list()

  if (missing(prior.beta)) prior.beta <- list()
  betadi <- bayessurvreg3.priorBeta(prior.beta, init, des)
  init$beta <- attr(betadi, "init")
  prior.beta <- attr(betadi, "prior.beta")
  
  if (missing(prior.beta2)) prior.beta2 <- list()
  betadi2 <- bayessurvreg3.priorBeta(prior.beta2, init2, des2)
  init2$beta <- attr(betadi2, "init")
  prior.beta2 <- attr(betadi2, "prior.beta")
  
  ## Priors and inits for rho (correlation coefficient between two random intercepts)
  if (missing(priorinit.Nb)){
    rho <- bayessurvreg3.checkrho(rho=rho, doubly=doubly)
    if (rho$type.update == "fixed.zero") version <- 3
    else                                 version <- 31
  }
  else{
    rho <- list(type.update = "fixed.zero", init=0, sigmaL=0.1)
    rho <- bayessurvreg3.checkrho(rho=rho, doubly=doubly)    
    version <- 32
  }  
  
  ## Priors and inits for random intercept
  if (version %in% c(3, 31)){
    if (missing(prior.b)) prior.b <- list()  
    reffdi <- bayessurvreg3.priorb(prior.b=prior.b, init=init, design=des, mcmc.par=mcmc.par)
    prior.b  <- attr(reffdi, "prior.b")
    init     <- attr(reffdi, "init")
    mcmc.par <- attr(reffdi, "mcmc.par")
  
    ## Priors and inits for random intercept
    if (missing(prior.b2)) prior.b2 <- list()  
    reffdi2 <- bayessurvreg3.priorb(prior.b=prior.b2, init=init2, design=des2, mcmc.par=mcmc.par2)
    prior.b2  <- attr(reffdi2, "prior.b")
    init2     <- attr(reffdi2, "init")
    mcmc.par2 <- attr(reffdi2, "mcmc.par")
  }
  else{
    if (version == 32){
      Reff <- bayessurvreg3.priorinitNb(priorinit.Nb=priorinit.Nb, init=init, init2=init2, design=des, design2=des2, doubly=doubly)
      reffdi  <- Reff$reffdi
      reffdi2 <- Reff$reffdi2      
      priorinit.Nb <- attr(Reff, "priorinit.Nb")
      init         <- attr(Reff, "init")
      init2        <- attr(Reff, "init2")
      rm(list="Reff")
    }
    else{
      stop("It is strange but I cannot determine the version to be used")
    }  
  }  
  
  ## Priors and inits for error G-spline, (censored) observations and observational allocations
  if (!doubly) prior2 <- list()
  prinit <- bayessurvreg3.priorInit(prior, init, des, mcmc.par, prior2, init2, des2, mcmc.par2, doubly)
  init      <- attr(prinit, "init")
  prior     <- attr(prinit, "prior")  
  mcmc.par  <- attr(prinit, "mcmc.par")
  init2     <- attr(prinit, "init2")
  prior2    <- attr(prinit, "prior2")  
  mcmc.par2 <- attr(prinit, "mcmc.par2")


  ## Compute quantities to determine the space needed to be allocated
  ##   and numbers of iterations in different phases
  if (nsimul$nburn >= nsimul$niter) nsimul$nburn <- nsimul$niter - 1
  if (nsimul$nburn < 0) nsimul$nburn <- 0
 
  if (nsimul$nburn == 0) nruns <- 1
  else                   nruns <- 2

  nrun <- numeric(2)
  nrun[2] <- nsimul$niter - nsimul$nburn
  nrun[1] <- nsimul$nburn

  nwrite.run <- nrun
  nwrite.run[nsimul$nwrite <= nrun] <- nsimul$nwrite   
  max.nwrite <- max(nwrite.run)

  ## Write headers to files with stored values  
  bayessurvreg3.writeHeaders(dir = dir, doubly = doubly, prior.init = prinit,
                             priorb.di = reffdi, priorb2.di = reffdi2, store = store, design = des, design2 = des2, version = version,
                             mclass = mclass)  

  ## Combine similar parameters into one vector  
  dims <- c(nobs, as.numeric(doubly))
  storeV <- c(store$a, store$y, store$r, store$b, store$a2, store$y2, store$r2, store$b2, store$a.b, store$r.b, store$a.b2, store$r.b2)
  nsimul.run1 <- c(nrun[1], nsimul$nthin, nwrite.run[1])
  nsimul.run2 <- c(nrun[2], nsimul$nthin, nwrite.run[2])
  names(nsimul.run1) <- names(nsimul.run2) <- c("niter", "nthin", "nwrite")
  nsample1 <- nsimul.run1["niter"] %/% nsimul.run1["nthin"]
  nsample2 <- nsimul.run2["niter"] %/% nsimul.run2["nthin"]  
  
  cat("Simulation started on                       ", date(), "\n", sep = "")
  fit <- .C("bayessurvreg2", as.character(dir),
                             dims = as.integer(dims),
                             X1 = as.double(if(des$nX) t(des$X) else 0),
                             X2 = as.double(if(des2$nX) t(des2$X) else 0),
                             y1.left = as.double(prinit$y.left),
                             y1.right = as.double(prinit$y.right),
                             status1 = as.integer(prinit$status),
                             t2.left = as.double(prinit$t2.left),
                             t2.right = as.double(prinit$t2.right),
                             status2 = as.integer(prinit$status2),
                             iPML = double(dims[1]),
                             Ys1 = as.double(prinit$y),
                             Ys2 = as.double(prinit$y2),
                             r1 = as.integer(prinit$r),
                             r2 = as.integer(prinit$r2),
                             specif = as.integer(prinit$specification),
                             r1.b = as.integer(reffdi$r),
                             r2.b = as.integer(reffdi2$r),
                             specif.b = as.integer(c(reffdi$specification, reffdi2$specification)),
                             GsplI1 = as.integer(prinit$Gparmi),
                             GsplD1 = as.double(prinit$Gparmd),
                             GsplI2 = as.integer(prinit$Gparmi2),
                             GsplD2 = as.double(prinit$Gparmd2),
                             priorBetaI1 = as.integer(betadi$parmI),
                             priorBetaD1 = as.double(betadi$parmD),
                             priorBetaI2 = as.integer(betadi2$parmI),
                             priorBetaD2 = as.double(betadi2$parmD),
                             priorbI1 = as.integer(reffdi$bparmI),
                             priorbD1 = as.double(reffdi$bparmD),
                             priorbI2 = as.integer(reffdi2$bparmI),
                             priorbD2 = as.double(reffdi2$bparmD),
                             priorCovMatI1 = as.integer(reffdi$GsplI),
                             priorCovMatD1 = as.double(reffdi$GsplD),
                             priorCovMatI2 = as.integer(reffdi2$GsplI),
                             priorCovMatD2 = as.double(reffdi2$GsplD),
                             rho             = as.double(rho$init),
                             rho.accept      = integer(nsample1),
                             rho.type.update = as.integer(rho$typeI),
                             rho.sigmaL      = as.double(rho$sigmaL),
                             mclass.sens.spec = as.double(c(mclass$init.sens, mclass$init.spec)),            
                             mclass.logtime   = as.double(mclass$logTime),
                             mclass.status    = as.integer(mclass$Status),
                             mclass.paramI    = as.integer(c(mclass$nModel, mclass$nExaminer, mclass$nFactor, mclass$nvisit, mclass$Examiner, mclass$Factor)),
                             mclass.paramD    = as.double(mclass$Prior),
                             iter = as.integer(prinit$iter),
                             nsimul = as.integer(nsimul.run1),
                             store = as.integer(storeV),
                             version = as.integer(version),
                             mainSimul = as.integer(0),            
                             err = integer(1),
           PACKAGE = thispackage)
  
  if (fit$err != 0) stop ("Something went wrong during the simulation.")
  cat("Burn-up finished on                         ", date(), "   (iteration ", fit$iter, ")", "\n", sep = "")
  
  ## Rewrite sampled values by new files
  bayessurvreg3.writeHeaders(dir = dir, doubly = doubly, prior.init = prinit,
                             priorb.di = reffdi, priorb2.di = reffdi2, store = store, design = des, design2 = des2, version = version,
                             mclass = mclass)  
  
  ## Main simulation
  fit <- .C("bayessurvreg2", as.character(dir),
                             dims = as.integer(dims),
                             X1 = as.double(fit$X1),
                             X2 = as.double(fit$X2),
                             y1.left = as.double(fit$y1.left),
                             y1.right = as.double(fit$y1.right),
                             status1 = as.integer(fit$status1),
                             t2.left = as.double(fit$t2.left),
                             t2.right = as.double(fit$t2.right),
                             status2 = as.integer(fit$status2),
                             iPML = double(dims[1]),            
                             Ys1 = as.double(fit$Ys1),
                             Ys2 = as.double(fit$Ys2),
                             r1 = as.integer(fit$r1),
                             r2 = as.integer(fit$r2),
                             specif = as.integer(fit$specif),
                             r1.b = as.integer(fit$r1.b),
                             r2.b = as.integer(fit$r2.b),
                             specif.b = as.integer(fit$specif.b),            
                             GsplI1 = as.integer(fit$GsplI1),
                             GsplD1 = as.double(fit$GsplD1),
                             GsplI2 = as.integer(fit$GsplI2),
                             GsplD2 = as.double(fit$GsplD2),
                             priorBetaI1 = as.integer(fit$priorBetaI1),
                             priorBetaD1 = as.double(fit$priorBetaD1),
                             priorBetaI2 = as.integer(fit$priorBetaI2),
                             priorBetaD2 = as.double(fit$priorBetaD2),
                             priorbI1 = as.integer(fit$priorbI1),
                             priorbD1 = as.double(fit$priorbD1),
                             priorbI2 = as.integer(fit$priorbI2),
                             priorbD2 = as.double(fit$priorbD2),
                             priorCovMatI1 = as.integer(fit$priorCovMatI1),
                             priorCovMatD1 = as.double(fit$priorCovMatD1),
                             priorCovMatI2 = as.integer(fit$priorCovMatI2),
                             priorCovMatD2 = as.double(fit$priorCovMatD2),
                             rho             = as.double(fit$rho),
                             rho.accept      = integer(nsample2),
                             rho.type.update = as.integer(fit$rho.type.update),
                             rho.sigmaL      = as.double(fit$rho.sigmaL),
                             mclass.sens.spec = as.double(fit$mclass.sens.spec),            
                             mclass.logtime   = as.double(mclass$logTime),
                             mclass.status    = as.integer(mclass$Status),
                             mclass.paramI    = as.integer(c(mclass$nModel, mclass$nExaminer, mclass$nFactor, mclass$nvisit, mclass$Examiner, mclass$Factor)),
                             mclass.paramD    = as.double(mclass$Prior),            
                             iter = as.integer(fit$iter),
                             nsimul = as.integer(nsimul.run2),
                             store = as.integer(fit$store),
                             version = as.integer(fit$version),
                             mainSimul = as.integer(1),
                             err = integer(1),
           PACKAGE = thispackage)
  
  if (fit$err != 0) stop ("Something went wrong during the simulation.")  
  cat("Simulation finished on                      ", date(), "   (iteration ", fit$iter, ")", "\n", sep = "")     
  
  toreturn <- fit$iter
  attr(toreturn, "call") <- call

  attr(toreturn, "init") <- init
  attr(toreturn, "prior") <- prior
  attr(toreturn, "prior.beta") <- prior.beta
  attr(toreturn, "mcmc.par") <- mcmc.par

  if (doubly){
    attr(toreturn, "init2") <- init2
    attr(toreturn, "prior2") <- prior2
    attr(toreturn, "prior.beta2") <- prior.beta2
    attr(toreturn, "mcmc.par2") <- mcmc.par2
  }

  if (version %in% c(3, 31)){
    attr(toreturn, "prior.b") <- prior.b
    if (doubly) attr(toreturn, "prior.b2") <- prior.b2
  }  
  if (version == 31){
    attr(toreturn, "rho.accept") <- fit$rho.accept
  }  

  if (version == 32){
    attr(toreturn, "priorinit.Nb") <- fit$priorinit.Nb
  }

  if (mclass$nModel){
    attr(toreturn, "prior.classification") <- mclass$Prior
    attr(toreturn, "init.sens") <- classParam$init.sens
    attr(toreturn, "init.spec") <- classParam$init.spec

    ##### Calculate model fit statistics
    ##### -----------------------------------
    
    ### Deviance sample
    devGJK <- scan(paste(dir, "/devianceGJK.sim", sep = ""), skip = 1, quiet = TRUE)
    D.bar <- mean(devGJK, na.rm = TRUE)

    ### Posterior means of sensitivities/specificities
    sensspec <- matrix(scan(paste(dir, "/sens_spec.sim", sep = ""), skip = 1, quiet = TRUE), ncol = 2 * mclass$nExaminer * mclass$nFactor, byrow = TRUE)
    colnames(sensspec) <- scan(paste(dir, "/sens_spec.sim", sep = ""), what = character(), nlines = 1, quiet = TRUE)
    sensspec.bar <- apply(sensspec, 2, mean, na.rm = TRUE)
    rm(list = "sensspec")

    ### Init for posterior means of linear predictors
    eta.bar <- rep(0, des$n)
    
    ### Posterior means of random effects
    if (des$nrandom > 0){
      if (store$b){        
        bb <- matrix(scan(paste(dir, "/b.sim", sep = ""), skip = 1, quiet = TRUE), ncol = des$ncluster, byrow = TRUE)
        bb.bar <- apply(bb, 2, mean, na.rm = TRUE)
        rm(list = "bb")
      }else{
        bb.bar <- rep(NA, des$ncluster)
        warning("to have DIC calculated, store$b must be TRUE.")
      }
      eta.bar <- eta.bar + rep(bb.bar, des$nwithin)
    }

    ### Posterior means of fixed effects
    if (des$nX){
      beta <- matrix(scan(paste(dir, "/beta.sim", sep = ""), skip = 1, quiet = TRUE), ncol = des$nX, byrow = TRUE)
      colnames(beta) <- scan(paste(dir, "/beta.sim", sep = ""), what = character(), nlines = 1, quiet = TRUE)
      beta.bar <- apply(beta, 2, mean, na.rm = TRUE)
      rm(list = "beta")
      eta.bar <- eta.bar + as.numeric(des$X %*% beta.bar)
    }  
    
    ### Posterior means of G-spline (error term) parameters
    gspline <- matrix(scan(paste(dir, "/gspline.sim", sep = ""), skip = 1, quiet = TRUE), ncol = 5, byrow = TRUE)
    colnames(gspline) <- scan(paste(dir, "/gspline.sim", sep = ""), what = character(), nlines = 1, quiet = TRUE)
    gspline.bar <- apply(gspline, 2, mean, na.rm = TRUE)
    gspline.sigmaSq.bar <- mean(gspline[, "sigma1"]^2, na.rm = TRUE)
    gspline.scaleSq.bar <- mean(gspline[, "scale1"]^2, na.rm = TRUE)    
    rm(list = "gspline")

    mweight <- matrix(scan(paste(dir, "/mweight.sim", sep = ""), skip = 1, quiet = TRUE), ncol = 2 * prinit$Gparmi["K1"] + 1, byrow = TRUE)
    mweight.bar <- apply(mweight, 2, mean, na.rm = TRUE)
    rm(list = "mweight")
    
    ### Input for D.in.bar to be calculated in C++
    if (any(is.na(eta.bar))){
      D.in.bar <- NA      
    }else{
      maxnvisit <- max(mclass$nvisit)
      nss <- mclass$nExaminer * mclass$nFactor
      inD <- .C("iPML_misclass_GJK",     
                iPML  = double(des$n),
                dwork = double(3 * (1 + maxnvisit)),
                min_etaM  = as.double((-1) * eta.bar),
                sens      = as.double(sensspec.bar[1:nss]),
                spec      = as.double(sensspec.bar[(nss + 1):(2 * nss)]),
                logvtime  = as.double(mclass$logTime),
                status    = as.integer(mclass$Status),
                nExaminer = as.integer(mclass$nExaminer),
                nFactor   = as.integer(mclass$nFactor),
                nvisit    = as.integer(mclass$nvisit),
                maxnvisit = as.integer(maxnvisit),
                Examiner  = as.integer(mclass$Examiner),
                Factor    = as.integer(mclass$Factor),
                gg_K      = as.integer(prinit$Gparmi["K1"]),
                gg_gamma  = as.double(gspline.bar["gamma1"]),
                gg_delta  = as.double(gspline.bar["delta1"]),
                gg_sigma  = as.double(sqrt(gspline.sigmaSq.bar)),
                #gg_sigma  = as.double(gspline.bar["sigma1"]),
                gg_intcpt = as.double(gspline.bar["intercept1"]),
                gg_scale  = as.double(sqrt(gspline.scaleSq.bar)),
                #gg_scale  = as.double(gspline.bar["scale1"]),
                gg_w      = as.double(mweight.bar),
                nP        = as.integer(des$n),
                PACKAGE = thispackage)
      D.in.bar <- -2 * sum(log(inD$iPML))      
    }
           
    ### DIC and log(PML)
    pD <- D.bar - D.in.bar
    DIC <- D.bar + pD
    logPML <- sum(log(fit$iPML))
    
    attr(toreturn, "fitStat") <- c(logPML, DIC, D.bar, D.in.bar, pD)
    names(attr(toreturn, "fitStat")) <- c("logPML", "DIC", "D.bar", "D.in.bar", "pD")
    attr(toreturn, "iPML")   <- fit$iPML
    
  }
  
  class(toreturn) <- "bayessurvreg3"
  return(toreturn)    
}
