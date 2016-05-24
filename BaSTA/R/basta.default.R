# FILE NAME:   basta.default
# AUTHOR:      Fernando Colchero
# DATE:        07/Nov/2015
# VERSION:     1.9.4
# DESCRIPTION: This function estimates age-specific mortality from capture-
#              recapture/recovery (CRR) data when a large proportion of (or all) 
#              the records have unknown times of birth and death. It uses the 
#              framework described by Colchero & Clark (2012) J Anim Ecol.
# COMMENTS:    None
# ============================================================================ #
basta.default <- function(object, studyStart, studyEnd, minAge = 0, 
    model = "GO", shape = "simple", covarsStruct = "fused", niter = 11000, 
    burnin = 1001, thinning = 20, recaptTrans = studyStart, 
    nsim = 1, parallel = FALSE, ncpus = 2, lifeTable = TRUE, 
    updateJumps = TRUE, ...) {
  
  argList <- list(...)
  bastaIntVars <- c("algObj", "defTheta", "CalcMort", "CalcSurv", 
      "dataObj", "covObj", "userPars", "fullParObj", "agesIni", 
      "parsIni", "priorAgeObj", "parsCovIni", "postIni", 
      "jumps", ".Random.seed")
  algObj <- .CreateAlgObj(model, shape, studyStart, studyEnd, minAge, 
      covarsStruct, recaptTrans, niter, burnin, thinning, updateJumps, nsim)
  .FindErrors(object, algObj)
  dataObj <- .PrepDataObj(object, algObj)
  defTheta <- .SetDefaultTheta(algObj)
  CalcMort <- .DefineMort(algObj)
  CalcSurv <- .DefineSurv(algObj)
  covObj <- .CreateCovObj(object, dataObj, algObj)
  algObj$covStruc <- class(covObj)[1]
  userPars <- .CreateUserPar(covObj, argList)
  fullParObj <- .BuildFullParObj(covObj, defTheta, algObj, userPars, 
      dataObj)
  agesIni <- .PrepAgeObj(dataObj, algObj)
  parsIni <- .DefineIniParObj(fullParObj)
  parsCovIni <- .BuildParCovObj(covObj, parsIni)
  priorAgeObj <- .SetPriorAgeDist(fullParObj, CalcSurv, dataObj, covObj, 
      parsIni, parsCovIni)
  postIni <- .BuildPostObj(agesIni, parsIni, parsCovIni, dataObj, 
      CalcSurv, priorAgeObj, fullParObj, covObj)
  jumps <- list()
  jumps$theta <- fullParObj$theta$jump
  if (fullParObj$class[1] == "theGam") {
    jumps$gamma <- fullParObj$gamma$jump
  }
  Start <- Sys.time()
  if(updateJumps) {
    .jumpObjIni <- .PrepJumpObj(fullParObj, covObj)
    updatedJumps <- .RunIniUpdJump(argList, algObj, defTheta, 
        CalcMort, CalcSurv, dataObj, covObj, userPars, fullParObj, 
        agesIni, parsIni, priorAgeObj, parsCovIni, postIni, jumps, 
        .jumpObjIni)
    jumps <- updatedJumps$updJumps
  } 
  if (nsim > 1) {
    cat("Multiple simulations started...\n\n") 
    if (parallel) {
			opp <- options()
			options(warn = -1)
			sfInit(parallel = TRUE, cpus = ncpus)
			sfExport(list = c(bastaIntVars, ".Random.seed"))
			sfLibrary("BaSTA", character.only = TRUE, 
					warn.conflicts = FALSE)
			bastaOut <- sfClusterApplyLB(1:nsim, .RunBastaMCMC, algObj,
					defTheta, CalcMort, CalcSurv, dataObj, covObj, userPars, 
					fullParObj, agesIni, parsIni, priorAgeObj, parsCovIni, postIni, 
					jumps)
			sfRemoveAll(hidden = TRUE)
			sfStop()
			options(opp)
		} else {
      bastaOut <- lapply(1:nsim, .RunBastaMCMC, algObj, defTheta, 
          CalcMort, CalcSurv, dataObj, covObj, userPars, fullParObj, 
          agesIni, parsIni, priorAgeObj, parsCovIni, postIni, jumps)
    }
  } else {
    cat("Simulation started...\n\n")
    bastaOut <- lapply(1:nsim, .RunBastaMCMC, algObj, defTheta, 
        CalcMort, CalcSurv, dataObj, covObj, userPars, fullParObj, 
        agesIni, parsIni, priorAgeObj, parsCovIni, postIni, jumps)
  }
  End <- Sys.time()
  names(bastaOut) <- paste("sim.", 1:nsim, sep = "")
  compTime <- round(as.numeric(End-Start, units = units(End - Start)), 2)
  cat(sprintf("Total MCMC computing time: %.2f %s.\n\n", compTime, 
          units(End - Start)))
  bastaResults <- .CalcDiagnost(bastaOut, algObj, covObj, defTheta, 
      fullParObj, dataObj, parsIni, parsCovIni, agesIni, postIni, CalcSurv, 
			priorAgeObj)
  bastaResults$settings <- c(niter, burnin, thinning, nsim)
  names(bastaResults$settings) <- c("niter", "burnin", "thinning", "nsim")
  bastaResults$modelSpecs <- 
    c(model, shape, covarsStruct, minAge,
        paste(names(covObj$cat), collapse = ", "), 
          paste(names(covObj$cont), collapse = ", "))
  names(bastaResults$modelSpecs) <- c("model", "shape", "Covar. structure", 
      "minAge", "Categorical", "Continuous")
  bastaResults <- .CalcQuants(bastaOut, bastaResults, defTheta, fullParObj, 
      algObj, dataObj, CalcSurv, CalcMort, covObj, agesIni)
  bastaResults$jumpPriors <- 
      cbind(c(jumps$theta, jumps$gamma), 
          c(fullParObj$theta$priorMean, fullParObj$gamma$priorMean), 
          c(fullParObj$theta$priorSd, fullParObj$gamma$priorSd))
  dimnames(bastaResults$jumpPriors) <- 
      list(fullParObj$allNames[substr(fullParObj$allNames, 1, 2) != "pi"],
          c("Jump.sds", "Prior.means", "Prior.sds"))
  bastaResults$parsForPlot <- list()
  for (pp in names(bastaOut)) {
    bastaResults$parsForPlot[[pp]] <- 
        bastaOut[[pp]]$par[seq(1, algObj$niter, algObj$thinning), ]
  }
  bastaResults$lifeTable <- .CalcLifeTable(bastaResults, lifeTable, object,
      covObj, algObj)
  bastaResults$version <- packageDescription("BaSTA")$Version
  # Define class for output object:
  class(bastaResults) <- "basta"
  return(bastaResults)
}

