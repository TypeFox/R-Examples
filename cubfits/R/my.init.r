### Get all specific functions according to the options.
### This is only called at the begining of each model.

my.init.function <- function(model = .CF.CT$model[1],
    # type.B = .CF.CT$type.B[1],
    type.p = .CF.CT$type.p[1],
    type.Phi = .CF.CT$type.Phi[1], model.Phi = .CF.CT$model.Phi[1],
    init.Phi = .CF.CT$init.Phi[1], init.fit = .CF.CT$init.fit[1],
    parallel = .CF.CT$parallel[1], adaptive = .CF.CT$adaptive[1]){
  ### Overwite conditions changed due to limited implementations.
  if(model.Phi == "logmixture"){
    type.p <- "logmixture"
  }
  if(.CF.CONF$estimate.bias.Phi){
    type.p <- "lognormal_bias"
  }
  if(type.p == "lognormal_bias"){
    .CF.CONF$estimate.bias.Phi <- TRUE
  }

  ### Check.
  if(any(parallel %in% c("task.pull", "pbdLapply"))){
    if(interactive()){
      stop("task.pull or pbdLapply is only for non-interactive.")
    }
    if(!is.loaded("spmd_initialize", PACKAGE = "pbdMPI")){
      stop("pbdMPI is not loaded.")
    }
  }

  ### Assign model.
  assign("model", model, envir = .cubfitsEnv)

  ### Coefficients.
  my.ncoef <- get.my.ncoef(model = model)
  my.coefnames <- get.my.coefnames(model = model)

  ### Major functions.
  my.drawBConditionalAll <- get.my.drawBConditionalAll(type = init.fit)
  my.pPropType <- get.my.pPropType(type = type.p)
  my.proposePhiAll <- get.my.proposePhiAll(type = type.Phi)
  my.fitMultinomOne <- get.my.fitMultinomOne(model = model)
  my.fitMultinomAll <- get.my.fitMultinomAll(parallel = parallel)
  my.logdmultinomCodOne <- get.my.logdmultinomCodOne(model = model)
  my.logdmultinomCodAllR <- get.my.logdmultinomCodAllR(parallel = parallel)
  my.logPosteriorAll <- get.my.logPosteriorAll(model.Phi = model.Phi)
  my.logLAll <- get.my.logLAll(model.Phi = model.Phi)

  ### Prediction functions.
  my.proposePhiAllPred <- get.my.proposePhiAllPred(type = type.Phi)
  my.logPosteriorAllPred <- get.my.logPosteriorAllPred(model.Phi = model.Phi)
  my.logLAllPred <- get.my.logLAllPred(model.Phi = model.Phi)
  my.estimatePhiOne <- get.my.estimatePhiOne(init.Phi = init.Phi, model = model)
  my.estimatePhiAll <- get.my.estimatePhiAll(parallel = parallel)

  ### No observation prior.
  my.pPropTypeNoObs <- get.my.pPropTypeNoObs(type = type.p)

  ### Utility functions.
  my.cat <- get.my.cat(parallel = parallel)
  my.print <- get.my.print(parallel = parallel)
  my.stop <- get.my.stop(parallel = parallel)
  my.dump <- get.my.dump(parallel = parallel)

  ### Acceptance rate and adaptive functions.
  my.update.DrawScale <- get.my.update.DrawScale(adaptive = adaptive)

  ### Return some information back to parent environment.
  ret <- list(my.ncoef = my.ncoef,
              ### for training.
              my.drawBConditionalAll = my.drawBConditionalAll,
              my.pPropType = my.pPropType,
              my.proposePhiAll =  my.proposePhiAll,
              my.fitMultinomOne = my.fitMultinomOne,
              my.fitMultinomAll = my.fitMultinomAll,
              my.logdmultinomCodOne = my.logdmultinomCodOne,
              my.logdmultinomCodAllR = my.logdmultinomCodAllR,
              my.logPosteriorAll = my.logPosteriorAll,
              my.logLAll = my.logLAll,
              ### for prediction.
              my.proposePhiAllPred =  my.proposePhiAllPred,
              my.logPosteriorAllPred = my.logPosteriorAllPred,
              my.logLAllPred = my.logLAllPred,
              my.estimatePhiOne = my.estimatePhiOne,
              my.estimatePhiAll = my.estimatePhiAll,
              ### for no observation.
              my.pPropTypeNoObs = my.pPropTypeNoObs,
              ### for utility.
              my.cat = my.cat,
              my.print = my.print,
              my.stop = my.stop,
              my.dump = my.dump,
              ### for adaptive function.
              my.update.DrawScale = my.update.DrawScale
             )
  invisible(ret)
} # End of my.init.function().

