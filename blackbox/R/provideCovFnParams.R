provideCovFnParams <- function(gcvKgpointls,
                             gcvnuniquerows,
                             fittedNames=blackbox.getOption("fittedNames"),
                             ycolname=blackbox.getOption("ycolname"),
                             minSmoothness=blackbox.getOption("minSmoothness"),
                             miscOptions=blackbox.getOption("miscOptions"),
                             initCovFnParam=NULL,
                             cleanResu=NULL,
                             verbosity=blackbox.getOption("verbosity"),
                             optimizers=blackbox.getOption("optimizers")
) {
  hglmLambdaEst <- hglmPhiEst <- lambdaEst <- NA
  if ("optimizeKriging" %innc% miscOptions) {
    Cfit <- CKrigcoefs(gcvKgpointls[, c(fittedNames, ycolname)], initCovFnParam=initCovFnParam,
                       nuniquerows=gcvnuniquerows,optimizers=optimizers)
    CovFnParam <- Cfit$covfnparam ##includes smoothness (FR->FR 09/2015: !! c'est le comportement de CKrigcoefs seulement quand option(minSmoothness) n'est pas nul...!!)
    if("HGLM" %innc% miscOptions) {
      hglmLambdaEst <- Cfit$hglmLambda
      hglmPhiEst <- Cfit$hglmPhi
    } else {
      lambdaEst <- Cfit$lambda
    }
    if (blackbox.getOption("verbosity")) {
      lllocalst <- paste("Estimation of covariance parameters required ", Cfit$fnEvalCount, " CV function evaluations")
      cat(lllocalst, "\n")
      lllocalst <- NA
    }
    names(CovFnParam) <- c(fittedNames, "smoothness") ##global
    if (blackbox.getOption("verbosity")) cat("Cross-validation estimates of correlation function parameters:", "\n")
    if( !is.null(cleanResu) && CovFnParam["smoothness"]<3.95) { ## FR->FR remplacer tous les tests null cleanResu par des tests class file...?
      if (minSmoothness>1.99) { ## ie if we expect high smoothness
        message.redef("(!!!) Estimated smoothness parameter <3.95: likelihood prediction may be very poor")
        write("(!) Estimated smoothness parameter <3.95: likelihood prediction may be very poor", file=cleanResu)
      }
    } ## FR->FR le mettre dans returnCode ?
  } else if (!(blackbox.getOption("CovFnParamInSettingsBool"))){ ## smoothness fixed in settings but covariance params not fixed in settings
    ## if we're here, then ( ! optimizeKriging) ... probably obsolete...
    cat("Correlation function parameters set by rough heuristics", "\n")
    CovFnParam <- ((blackbox.getOption("FONKgUp")-blackbox.getOption("FONKgLow"))[fittedNames])/blackbox.getOption("metarange")
    CovFnParam <- c(CovFnParam, unlist(list(smoothness=minSmoothness)))
  } else { ## everything is in settings
    cat("Correlation function parameters set by user:", "\n")
    CovFnParam <- c(blackbox.getOption("CovFnParam"), unlist(list(smoothness=minSmoothness)))
  }
  names(CovFnParam) <- c(fittedNames, "smoothness")
  locNvalues <- CovFnParam
  userNames <- sapply(fittedNames, formatName, format="ASCII")
  names(locNvalues) <- c(userNames, "smoothness") ##global
  if (blackbox.getOption("verbosity")) {
    print(locNvalues)
    cat("\n")
  }
  nuggetNA <- ( (("HGLM" %innc% miscOptions) && is.na(hglmLambdaEst))
                || ( ( ! ("HGLM" %innc% miscOptions)) && is.na(lambdaEst)))
  if (nuggetNA) {
    ## CovFnParam In Settings => CKrigcoefs has not yet been run and lambda is still NA
    Cfit <- CKrigcoefs(gcvKgpointls[, c(fittedNames, ycolname)],
                       nuniquerows=gcvnuniquerows,
                       covfnparamA=CovFnParam,
                       lambdaA=NA,
                       optimizers=optimizers)
    if ("HGLM" %innc% miscOptions) {
      hglmLambdaEst <- Cfit$hglmLambda
      hglmPhiEst <- Cfit$hglmPhi
    } else {
      lambdaEst <- Cfit$lambda
    }
  }
  return(list(CovFnParam=CovFnParam, lambdaEst=lambdaEst,method=Cfit$method))
}
