#' @rdname weights-methods
#' @rdname fmrsfit-class
#' @aliases weights,weights-method
setMethod("weights", signature = "fmrsfit", weights.fmrsfit)

#' @rdname residuals-methods
#' @rdname fmrsfit-class
#' @aliases residuals,residuals-method
setMethod("residuals", signature = "fmrsfit", residuals.fmrsfit)

#' @rdname nobs-methods
#' @rdname fmrsfit-class
#' @aliases nobs,nobs-method
setMethod("nobs", signature = "fmrsfit", nobs.fmrsfit)

#' @rdname ncov-methods
#' @rdname fmrsfit-class
#' @aliases ncov,ncov-method
setMethod("ncov", signature = "fmrsfit", ncov.fmrsfit)

#' @rdname ncomp-methods
#' @rdname fmrsfit-class
#' @aliases ncomp,ncomp-method
setMethod("ncomp", signature = "fmrsfit", ncomp.fmrsfit)

#' @rdname mixProp-methods
#' @rdname fmrsfit-class
#' @aliases mixProp,mixProp-method
setMethod("mixProp", signature = "fmrsfit", mixProp.fmrsfit)

#' @rdname logLik-methods
#' @rdname fmrsfit-class
#' @aliases logLik,logLik-method
setMethod("logLik", signature = "fmrsfit", logLik.fmrsfit)

#' @rdname fitted-methods
#' @rdname fmrsfit-class
#' @aliases fitted,fitted-method
setMethod("fitted", signature = "fmrsfit", fitted.fmrsfit)

#' @rdname dispersion-methods
#' @rdname fmrsfit-class
#' @aliases dispersion,dispersion-method
setMethod("dispersion", signature = "fmrsfit", dispersion.fmrsfit)

#' @rdname coefficients-methods
#' @rdname fmrsfit-class
#' @aliases coefficients,coefficients-method
setMethod("coefficients", signature = "fmrsfit", coefficients.fmrsfit)

#' @rdname BIC-methods
#' @rdname fmrsfit-class
#' @aliases BIC,BIC-method
setMethod("BIC", signature = "fmrsfit", BIC.fmrsfit)

#' @rdname summary-methods
#' @rdname fmrsfit-class
#' @aliases summary,summary-method
setMethod("summary", signature = "fmrsfit", summary.fmrsfit)

#' @rdname show-methods
#' @rdname fmrsfit-class
#' @aliases show,show-method
setMethod("show", signature = "fmrsfit", show.fmrsfit)


#' @rdname fmrs.mle-methods
#' @aliases fmrs.mle-method
setMethod(f="fmrs.mle", definition=function(y,
                                         delta,
                                         x,
                                         nComp = 2,
                                         disFamily = "lnorm",
                                         initCoeff,
                                         initDispersion,
                                         initmixProp,
                                         lambRidge = 0,
                                         nIterEM = 400,
                                         nIterNR = 2,
                                         conveps = 1e-8,
                                         convepsEM = 1e-8,
                                         convepsNR = 1e-8,
                                         porNR = 2){
  if(missing(y) | !is.numeric(y))
    stop("A numeric response vector must be provided.")
  if(missing(x) | !is.numeric(x))
    stop("A numeric matrix for covariates must be provided.")
  if(missing(delta) & (disFamily!="norm"))
    stop("A censoring indicator vector with 0 or 1 values must be provided.")

    if((nComp<2) ) {
    stop("An interger greater than 2 for the order of mixture model
         must be provided.")
  }
  nCov = dim(x)[2]
  n = length(y)
  if(missing(initCoeff)) initCoeff = rnorm((nCov+1)*nComp)
  if(missing(initDispersion)) initDispersion = rep(1,nComp)
  if(missing(initmixProp)) initmixProp = rep(1/nComp,nComp)

  coef0 <- matrix(c(initCoeff), nrow = nComp, ncol = nCov+1, byrow = TRUE)

  if(is.null(colnames(x))){
    xnames <- c("Intercept",c(paste("Cov",1:nCov,sep=".")))
  } else{
    xnames <- c("Intercept",colnames(x))
  }
  comnames <- c(paste("Comp",1:nComp,sep="."))

  if(disFamily == "norm"){
    model = "FMR"
    delta = rep(1, n)
    res=.C("FMR_Norm_Surv_EM_MLE", PACKAGE="fmrs",
           y = as.double(y),
           x = as.double(as.vector(unlist(x))),
           delta = as.double(delta),
           Lambda.Ridge = as.double(lambRidge),
           Num.Comp = as.integer(nComp),
           Num.Cov = as.integer(nCov),
           Sample.Size = as.integer(n),
           Num.iteration = as.integer(nIterEM),
           Max.iterEM.used = as.integer(0),
           Initial.Intercept = as.double(unlist(coef0[,1])),
           Initial.Coefficient = as.double(unlist(t(coef0[,-1]))),
           Initial.Dispersion = as.double(initDispersion),
           Initial.mixProp = as.double(initmixProp),
           conv.eps = as.double(conveps),
           conv.eps.em = as.double(convepsEM),
           Intecept.Hat = as.double(rep(0,nComp)),
           Coefficient.Hat = as.double(rep(0,nComp*nCov)),
           Dispersion.Hat = as.double(rep(0,nComp)),
           mixProp.Hat = as.double(rep(0,nComp)),
           LogLikelihood = as.double(0),
           BIC = as.double(0),
           AIC = as.double(0),
           GCV = as.double(0),
           EBIC1 = as.double(0),
           EBIC5 = as.double(0),
           GIC = as.double(0),
           predict = as.double(rep(0,n*nComp)),
           residual = as.double(rep(0,n*nComp)),
           tau = as.double(rep(0,n*nComp))
    )
  }else if(disFamily == "lnorm"){
    model = "FMAFTR"
    logy = log(y)
    res=.C("FMR_Norm_Surv_EM_MLE", PACKAGE="fmrs",
           y = as.double(logy),
           x = as.double(as.vector(unlist(x))),
           delta = as.double(delta),
           Lambda.Ridge = as.double(lambRidge),
           Num.Comp = as.integer(nComp),
           Num.Cov = as.integer(nCov),
           Sample.Size = as.integer(n),
           Num.iteration = as.integer(nIterEM),
           Max.iterEM.used = as.integer(0),
           Initial.Intercept = as.double(unlist(coef0[,1])),
           Initial.Coefficient = as.double(unlist(t(coef0[,-1]))),
           Initial.Dispersion = as.double(initDispersion),
           Initial.mixProp = as.double(initmixProp),
           conv.eps = as.double(conveps),
           conv.eps.em = as.double(convepsEM),
           Intecept.Hat = as.double(rep(0,nComp)),
           Coefficient.Hat = as.double(rep(0,nComp*nCov)),
           Dispersion.Hat = as.double(rep(0,nComp)),
           mixProp.Hat = as.double(rep(0,nComp)),
           LogLikelihood = as.double(0),
           BIC = as.double(0),
           AIC = as.double(0),
           GCV = as.double(0),
           EBIC1 = as.double(0),
           EBIC5 = as.double(0),
           GIC = as.double(0),
           predict = as.double(rep(0,n*nComp)),
           residual = as.double(rep(0,n*nComp)),
           tau = as.double(rep(0,n*nComp))
    )

  }else if(disFamily == "weibull"){
    model = "FMAFTR"
    logy = log(y)
    res=.C("FMR_Weibl_Surv_EM_MLE", PACKAGE="fmrs",
           y = as.double(logy),
           x = as.double(as.vector(unlist(x))),
           delta = as.double(delta),
           Lambda.Ridge = as.double(lambRidge),
           Num.Comp = as.integer(nComp),
           Num.Cov = as.integer(nCov),
           Sample.Size = as.integer(n),
           Num.iterationEM = as.integer(nIterEM),
           Num.iterationNR = as.integer(nIterNR),
           PortionNF = as.integer(porNR),
           Max.iterEM.used = as.integer(0),
           Initial.Intercept = as.double(c(coef0[,1])),
           Initial.Coefficient = as.double(c(t(coef0[,-1]))),
           Initial.Dispersion = as.double(initDispersion),
           Initial.mixProp = as.double(initmixProp),
           conv.eps.em = as.double(convepsEM),
           Intecept.Hat = as.double(rep(0,nComp)),
           Coefficient.Hat = as.double(rep(0,nComp*nCov)),
           Dispersion.Hat = as.double(rep(0,nComp)),
           mixProp.Hat = as.double(rep(0,nComp)),
           LogLikelihood = as.double(0),
           BIC = as.double(0),
           AIC = as.double(0),
           GCV = as.double(0),
           EBIC1 = as.double(0),
           EBIC5 = as.double(0),
           GIC = as.double(0),
           predict = as.double(rep(0,n*nComp)),
           residual = as.double(rep(0,n*nComp)),
           tau = as.double(rep(0,n*nComp))
    )
  }else{
    stop("The family of sub-distributions is not specified correctly.")
  }

  fit <- new("fmrsfit",
             y = y,
             delta = delta,
             x = x,
             nobs = n,
             ncov = nCov,
             ncomp = nComp,
             coefficients = array(rbind(res$Intecept.Hat,
                                        matrix(res$Coefficient.Hat,
                                               nrow = nCov, byrow = FALSE)),
                                  dim = c(nCov+1, nComp),
                                  dimnames = list(xnames,comnames)),
             dispersion = array(res$Dispersion.Hat, dim =
                                c(1,nComp),dimnames = list(NULL,comnames)),
             mixProp = array(res$mixProp.Hat, dim =
                               c(1,nComp),dimnames = list(NULL,comnames)),
             logLik = res$LogLikelihood,
             BIC = res$BIC,
             nIterEMconv = res$Max.iterEM.used,
             disFamily = disFamily,
             lambRidge = lambRidge,
             model = model,
             fitted = array(matrix(res$predict, nrow = n, byrow = FALSE),
                            dim = c(n, nComp), dimnames =
                              list(NULL,comnames)),
             residuals = array(matrix(res$residual, nrow =n, byrow = FALSE),
                               dim = c(n, nComp),
                               dimnames = list(NULL,comnames)),
             weights = array(matrix(res$tau, nrow = n, byrow = FALSE),
                             dim = c(n, nComp), dimnames =
                               list(NULL,comnames))
  )
  return(fit)
})




#' @rdname fmrs.tunsel-methods
#' @aliases fmrs.tunsel-method
setMethod(f="fmrs.tunsel", definition=function(y,
                                               delta,
                                               x,
                                               nComp,
                                               disFamily = "lnorm",
                                               initCoeff,
                                               initDispersion,
                                               initmixProp,
                                               penFamily = "lasso",
                                               lambRidge = 0,
                                               nIterEM = 2000,
                                               nIterNR = 2,
                                               conveps = 1e-8,
                                               convepsEM = 1e-8,
                                               convepsNR = 1e-8,
                                               porNR = 2,
                                               gamMixPor = 1){

  if(missing(y) | !is.numeric(y))
    stop("A numeric response vector must be provided.")
  if(missing(x) | !is.numeric(x))
    stop("A numeric matrix for covariates must be provided.")
  if(missing(delta) & (disFamily!="norm"))
    stop("A censoring indicator vector with 0 or 1 values must be provided.")

  if((nComp<2) ) {
    stop("An interger greater than 2 for the order of mixture model
         must be provided.")
  }
  nCov = dim(x)[2]
  n = length(y)
  if(missing(initCoeff)) initCoeff = rnorm((nCov+1)*nComp)
  if(missing(initDispersion)) initDispersion = rep(1,nComp)
  if(missing(initmixProp)) initmixProp = rep(1/nComp,nComp)

  coef0 <- matrix(c(initCoeff), nrow = nComp, ncol = nCov+1, byrow = TRUE)

  if(penFamily == "lasso") myPenaltyFamily = 1
  else if (penFamily == "scad") myPenaltyFamily = 2
  else if (penFamily == "mcp") myPenaltyFamily = 3
  else if (penFamily == "sica") myPenaltyFamily = 4
  else if (penFamily == "adplasso") myPenaltyFamily = 5
  else if (penFamily == "hard") myPenaltyFamily = 6
  else {
    stop("Penalty is not correctly specified.")
    }

  if(disFamily == "norm"){
    model = "FMR"
    delta = rep(1, n)

    res=.C("FMR_Norm_Surv_CwTuneParSel", PACKAGE="fmrs",
           y = as.double(y),
           x = as.double(as.vector(unlist(x))),
           delta = as.double(delta),
           myPenaltyFamily = as.integer(myPenaltyFamily),
           Lambda.Ridge = as.double(lambRidge),
           Num.Comp = as.integer(nComp),
           Num.Cov = as.integer(nCov),
           Sample.Size = as.integer(n),
           Initial.Intercept = as.double(c(coef0[,1])),
           Initial.Coefficient = as.double(c(t(coef0[,-1]))),
           Initial.Dispersion = as.double(initDispersion),
           Initial.mixProp = as.double(initmixProp),
           conv.eps = as.double(conveps),
           conv.eps.em = as.double(convepsEM),
           GamMixPortion = as.double(gamMixPor),
           Opt.Lambda = as.double(rep(0,nComp))
    )
  }else if(disFamily == "lnorm"){
    model = "FMAFTR"
    logy = log(y)
    res=.C("FMR_Norm_Surv_CwTuneParSel", PACKAGE="fmrs",
           y = as.double(logy),
           x = as.double(as.vector(unlist(x))),
           delta = as.double(delta),
           myPenaltyFamily = as.integer(myPenaltyFamily),
           Lambda.Ridge = as.double(lambRidge),
           Num.Comp = as.integer(nComp),
           Num.Cov = as.integer(nCov),
           Sample.Size = as.integer(n),
           Initial.Intercept = as.double(c(coef0[,1])),
           Initial.Coefficient = as.double(c(t(coef0[,-1]))),
           Initial.Dispersion = as.double(initDispersion),
           Initial.mixProp = as.double(initmixProp),
           conv.eps = as.double(conveps),
           conv.eps.em = as.double(convepsEM),
           GamMixPortion = as.double(gamMixPor),
           Opt.Lambda = as.double(rep(0,nComp))
    )
  }else if(disFamily == "weibull"){
    model = "FMAFTR"
    logy = log(y)

    res=.C("FMR_Weibl_Surv_CwTuneParSel", PACKAGE = "fmrs",
           y = as.double(logy),
           x = as.double(as.vector(unlist(x))),
           delta = as.double(delta),
           myPenaltyFamily = as.integer(myPenaltyFamily),
           Lambda.Ridge = as.double(lambRidge),
           Num.Comp = as.integer(nComp),
           Num.Cov = as.integer(nCov),
           Sample.Size = as.integer(n),
           Initial.Intercept = as.double(c(coef0[,1])),
           Initial.Coefficient = as.double(c(t(coef0[,-1]))),
           Initial.Dispersion = as.double(initDispersion),
           Initial.mixProp = as.double(initmixProp),
           Num.NRiteration = as.double(nIterNR),
           Num.PortionNF = as.double(porNR),
           conv.eps = as.double(convepsNR),
           GamMixPortion = as.double(gamMixPor),
           Opt.Lambda = as.double(rep(0,nComp))
    )
  }else{
    stop("The family of sub-distributions is not specified correctly.")
  }


  lambdafit <- new("fmrstunpar",
                   ncomp = nComp,
                   lambPen = array(res$Opt.Lambda, dim = c(1,nComp),
                                   dimnames = c(list(NULL,
                                                     c(paste("Comp", 1:nComp,
                                                             sep = "."))))),
                   lambRidge = lambRidge,
                   disFamily = disFamily,
                   penFamily = penFamily,
                   model = model
  )
  return(lambdafit)
}
)



#' @rdname fmrs.varsel-methods
#' @aliases fmrs.varsel-method
setMethod(f="fmrs.varsel", definition=function(y,
                                               delta,
                                               x,
                                               nComp,
                                               disFamily = "lnorm",
                                               initCoeff,
                                               initDispersion,
                                               initmixProp,
                                               penFamily = "lasso",
                                               lambPen,
                                               lambRidge = 0,
                                               nIterEM = 2000,
                                               nIterNR = 2,
                                               conveps = 1e-8,
                                               convepsEM = 1e-8,
                                               convepsNR = 1e-8,
                                               porNR = 2,
                                               gamMixPor = 1){
  if(missing(y) | !is.numeric(y))
    stop("A numeric response vector must be provided.")
  if(missing(x) | !is.numeric(x))
    stop("A numeric matrix for covariates must be provided.")
  if(missing(delta) & (disFamily!="norm"))
    stop("A censoring indicator vector with 0 or 1 values must be provided.")

  if((nComp<2) ) {
    stop("An interger greater than 2 for the order of mixture model
         must be provided.")
  }
  nCov = dim(x)[2]
  n = length(y)
  if(missing(initCoeff)) initCoeff = rnorm((nCov+1)*nComp)
  if(missing(initDispersion)) initDispersion = rep(1,nComp)
  if(missing(initmixProp)) initmixProp = rep(1/nComp,nComp)

  coef0 <- matrix(c(initCoeff), nrow = nComp, ncol = nCov+1, byrow = TRUE)

  if(is.null(colnames(x))){
    xnames <- c("Intercept",c(paste("Cov",1:nCov,sep=".")))
  } else{
    xnames <- c("Intercept",colnames(x))
  }
  comnames <- c(paste("Comp",1:nComp,sep="."))


  if(penFamily == "lasso") myPenaltyFamily = 1
  else if (penFamily == "scad") myPenaltyFamily = 2
  else if (penFamily == "mcp") myPenaltyFamily = 3
  else if (penFamily == "sica") myPenaltyFamily = 4
  else if (penFamily == "adplasso") myPenaltyFamily = 5
  else if (penFamily == "hard") myPenaltyFamily = 6
  else {stop("Penalty is not correctly specified.")  }

  if(disFamily == "norm"){
    model = "FMR"
    delta = rep(1, n)

    res=.C("FMR_Norm_Surv_EM_VarSel", PACKAGE="fmrs",
           y = as.double(y),
           x = as.double(as.vector(unlist(x))),
           delta = as.double(delta),
           myPenaltyFamily = as.integer(myPenaltyFamily),
           Lambda.Pen = as.double(lambPen),
           Lambda.Ridge = as.double(lambRidge),
           Num.Comp = as.integer(nComp),
           Num.Cov = as.integer(nCov),
           Sample.Size = as.integer(n),
           NumIterationEM = as.integer(nIterEM),
           Max.iterEM.used = as.integer(0),
           Initial.Intercept = as.double(c(coef0[,1])),
           Initial.Coefficient = as.double(c(t(coef0[,-1]))),
           Initial.Dispersion = as.double(initDispersion),
           Initial.mixProp = as.double(initmixProp),
           conv.eps = as.double(conveps),
           conv.eps.em = as.double(convepsEM),
           GamMixPortion = as.double(gamMixPor),
           Intecept.Hat = as.double(rep(0,nComp)),
           Coefficient.Hat = as.double(rep(0,nComp*nCov)),
           Dispersion.Hat = as.double(rep(0,nComp)),
           mixProp.Hat = as.double(rep(0,nComp)),
           LogLikelihood = as.double(0),
           BIC = as.double(0),
           AIC = as.double(0),
           GCV = as.double(0),
           EBIC1 = as.double(0),
           EBIC5 = as.double(0),
           GIC = as.double(0),
           predict = as.double(rep(0,n*nComp)),
           residual = as.double(rep(0,n*nComp)),
           tau = as.double(rep(0,n*nComp))
    )

  }else if(disFamily == "lnorm"){
    model = "FMAFTR"
    logy = log(y)

    res=.C("FMR_Norm_Surv_EM_VarSel", PACKAGE="fmrs",
           y = as.double(logy),
           x = as.double(as.vector(unlist(x))),
           delta = as.double(delta),
           myPenaltyFamily = as.integer(myPenaltyFamily),
           Lambda.Pen = as.double(lambPen),
           Lambda.Ridge = as.double(lambRidge),
           Num.Comp = as.integer(nComp),
           Num.Cov = as.integer(nCov),
           Sample.Size = as.integer(n),
           NumIterationEM = as.integer(nIterEM),
           Max.iterEM.used = as.integer(0),
           Initial.Intercept = as.double(c(coef0[,1])),
           Initial.Coefficient = as.double(c(t(coef0[,-1]))),
           Initial.Dispersion = as.double(initDispersion),
           Initial.mixProp = as.double(initmixProp),
           conv.eps = as.double(conveps),
           conv.eps.em = as.double(convepsEM),
           GamMixPortion = as.double(gamMixPor),
           Intecept.Hat = as.double(rep(0,nComp)),
           Coefficient.Hat = as.double(rep(0,nComp*nCov)),
           Dispersion.Hat = as.double(rep(0,nComp)),
           mixProp.Hat = as.double(rep(0,nComp)),
           LogLikelihood = as.double(0),
           BIC = as.double(0),
           AIC = as.double(0),
           GCV = as.double(0),
           EBIC1 = as.double(0),
           EBIC5 = as.double(0),
           GIC = as.double(0),
           predict = as.double(rep(0,n*nComp)),
           residual = as.double(rep(0,n*nComp)),
           tau = as.double(rep(0,n*nComp))
    )
  }else if(disFamily == "weibull"){
    model = "FMAFTR"
    logy = log(y)

    res=.C("FMR_Weibl_Surv_EM_VarSel", PACKAGE="fmrs",
           y = as.double(logy),
           x = as.double(as.vector(unlist(x))),
           delta = as.double(delta),
           myPenaltyFamily = as.integer(myPenaltyFamily),
           Lambda.Pen = as.double(lambPen),
           Lambda.Ridge = as.double(lambRidge),
           Num.Comp = as.integer(nComp),
           Num.Cov = as.integer(nCov),
           Sample.Size = as.integer(n),
           NumIterationEM = as.integer(nIterEM),
           NumIterationNR = as.integer(nIterNR),
           PortionNF = as.integer(porNR),
           Max.iterEM.used = as.integer(0),
           Initial.Intercept = as.double(c(coef0[,1])),
           Initial.Coefficient = as.double(c(t(coef0[,-1]))),
           Initial.Dispersion = as.double(initDispersion),
           Initial.mixProp = as.double(initmixProp),
           conv.eps = as.double(convepsNR),
           conv.eps.em = as.double(convepsEM),
           GamMixPortion = as.double(gamMixPor),
           Intecept.Hat = as.double(rep(0,nComp)),
           Coefficient.Hat = as.double(rep(0,nComp*nCov)),
           Dispersion.Hat = as.double(rep(0,nComp)),
           mixProp.Hat = as.double(rep(0,nComp)),
           LogLikelihood = as.double(0),
           BIC = as.double(0),
           AIC = as.double(0),
           GCV = as.double(0),
           EBIC1 = as.double(0),
           EBIC5 = as.double(0),
           GIC = as.double(0),
           predict = as.double(rep(0,n*nComp)),
           residual = as.double(rep(0,n*nComp)),
           tau = as.double(rep(0,n*nComp))
    )
  }else{
    stop("The family of sub-distributions is not specified correctly.")
  }

  fit <- new("fmrsfit", y = y,
             delta = delta,
             x = x,
             nobs = n,
             ncov = nCov,
             ncomp = nComp,
             coefficients = array(rbind(res$Intecept.Hat,
                                        matrix(res$Coefficient.Hat,
                                               nrow = nCov, byrow = FALSE)),
                                  dim = c(nCov+1, nComp),
                                  dimnames = list(xnames,comnames)),
             dispersion = array(res$Dispersion.Hat, dim = c(1,nComp),dimnames =
                                list(NULL,comnames)),
             mixProp = array(res$mixProp.Hat, dim = c(1,nComp),dimnames =
                               list(NULL,comnames)),
             logLik = res$LogLikelihood,
             BIC = res$BIC,
             nIterEMconv = res$Max.iterEM.used,
             disFamily = disFamily,
             penFamily = penFamily,
             lambPen = array(lambPen, dim = c(1,nComp),dimnames =
                               list(NULL,comnames)),
             model = model,
             fitted = array(matrix(res$predict, nrow = n, byrow = FALSE),
                            dim = c(n, nComp), dimnames =
                              list(NULL,comnames)),
             residuals = array(matrix(res$residual, nrow =n, byrow = FALSE),
                               dim = c(n, nComp), dimnames =
                                 list(NULL,comnames)),
             weights = array(matrix(res$tau, nrow = n, byrow = FALSE),
                             dim = c(n, nComp), dimnames =
                               list(NULL,comnames))
  )
  return(fit)
}

)





#' @rdname fmrs.gendata-methods
#' @aliases fmrs.gendata-method
setMethod(f="fmrs.gendata", definition=function(nObs,
                                                nComp,
                                                nCov,
                                                coeff,
                                                dispersion,
                                                mixProp,
                                                rho,
                                                umax,
                                                disFamily = "lnorm")

{
  if(missing(disFamily)) disFamily = "lnorm"

  if(sum(mixProp) != 1)
    stop("The sum of mixing proportions must be 1.")
  if(sum(dispersion <= 0) != 0)
    stop("Dispersion parameters cannot be zero or negative.")
  if(rho > 1 | rho < -1)
    stop("The correlation cannot be less than -1 or greater thatn 1.")

  mu <- rep(0, nCov)
  Sigma <- diag(nCov)

  for(i in 1:nCov){
    for(j in 1:nCov){
      Sigma[i,j] <- rho^abs(i-j)
    }}

  X <- matrix(rnorm(nCov * nObs), nObs)
  X <- scale(X, TRUE, FALSE)
  X <- X %*% svd(X, nu = 0)$v
  X <- scale(X, FALSE, TRUE)
  eS <- eigen(Sigma, symmetric = TRUE)
  ev <- eS$values
  X <- drop(mu) + eS$vectors %*% diag(sqrt(pmax(ev, 0)), nCov) %*% t(X)
  nm <- names(mu)
  if (is.null(nm) && !is.null(dn <- dimnames(Sigma)))
    nm <- dn[[1L]]
  dimnames(X) <- list(nm, NULL)
  if (nObs == 1)
    cX = drop(X)
  else cX = t(X)
  cX <- scale(cX)
  colnames(cX) <-  paste("X", 1:nCov,sep = ".")

  coef0 <- matrix(coeff, nrow = nComp, ncol = nCov+1, byrow = TRUE)
  mixProp0 <- cumsum(mixProp)

  yobs <-c()
  c <- rep()
  dlt <- c()
  u <- c()
  tobs <- c()

  if(disFamily == "lnorm"){
    for(i in 1:nObs){
      epss <- rnorm(1)
      u1 <- runif(1)
      k = length(which(mixProp0<=u1)) + 1
      u[i] = k
      yobs[i] <- coef0[k,1] + coef0[k,-1] %*% cX[i,] + dispersion[k] * epss

      c[i] <- log(runif(1, 0, umax))
      tobs[i] <- exp(min(yobs[i],c[i]))
      dlt[i] <- (yobs[i] < c[i])*1
    }
  }else if(disFamily=="norm"){
    for(i in 1:nObs){
      epss <- rnorm(1)
      u1 <- runif(1)
      k = length(which(mixProp0<=u1)) + 1
      u[i] = k
      yobs[i] <- coef0[k,1] + coef0[k,-1] %*% cX[i,] + dispersion[k] * epss
      tobs[i] <- yobs[i]
      dlt[i] <- 1
    }
  }else if(disFamily=="weibull"){
    for(i in 1:nObs){
      ext <- log(rexp(1))
      u1 <- runif(1)
      k = length(which(mixProp0<=u1)) + 1
      yobs[i] <- coef0[k,1] + coef0[k,-1] %*% cX[i,] + dispersion[k] * ext

      c[i]<- log(runif(1, 0, umax))
      tobs[i] <- exp(min(yobs[i],c[i]))
      dlt[i] <- (yobs[i] < c[i])*1
    }
  } else{
    stop("The family of sub-distributions are not specified correctly.")
  }
  return(list(y = tobs, delta = dlt, x = cX, disFamily = disFamily))
}
)
