# This function estimates any of the the NA parameters, or returns c and d if all parameters are given
CKrigcoefs <- function(xy,
                     nuniquerows, ## as given by selectFn
                     covfnparamA=NA,
                     lambdaA=NA,
                     minSmoothness=blackbox.getOption("minSmoothness"),
                     maxSmoothness=blackbox.getOption("maxSmoothness"),
                     miscOptions=blackbox.getOption("miscOptions"),
                     initCovFnParam=NULL,
                     verbosity=blackbox.getOption("verbosity"),
                     optimizers=blackbox.getOption("optimizers")
){
  nrowxy <- nrow(xy)
  ncolxy <- ncol(xy)
  optimise <- any(is.na(c(covfnparamA,lambdaA)))
  method <- intersect(optimizers,c("bobyqa","L-BFGS-B","lbfgsb3")) ## non-default methods
  if(length(method)==0L) method <- "NLOPT_LN_BOBYQA" ## default for this part of code.
  if(length(method)!=1L) stop("length(method)!=1L in CKrigcoefs") ## incompatible with switch below
  verbosityFromObjective <- as.integer((optimise && "L-BFGS-B" %in% method))*verbosity ## FR->FR should control in the C++ code using the eval counter
  verbosityFromoptimizer <- switch(method,
                                   "bobyqa"=round(2*as.integer(optimise)*verbosity), ## default=2
                                   "L-BFGS-B"=0, ## verbosity from optim() is awkward
                                   as.integer(optimise)*verbosity
  )
  success <- newCSmooth(xy= t(xy), ## newCSmooth is an Rcpp export
             nrowxy=nrowxy,
             ncolxy=ncolxy,
             nuniquerows=nuniquerows,
             GCV=0, # 0 for GCV, 1 for match_fs2hat_pure_error
             optimiseBool=optimise,
             verbosity=verbosityFromObjective)
  if ( ! success ) {
    resu <- list()
    ## go directly  to deleteCSmooth
  } else if ( ! optimise ) {
    resu <- Krig_coef_Wrapper(covfnparamA,lambdaA)[c("c","d","CKrigidx")] # throwing away u and D
    # do not deleteCSmooth(), CSmooth pointer stored in CKrigptrTable
  } else {
    if (any(is.na(covfnparamA))) { ## estim covfnparamA
      if ( is.null(initCovFnParam) ) initCovFnParam <- rep(0, ncolxy)
      minSmoothness <- min(maxSmoothness,max(minSmoothness,1.001)); # >1, cf def Matern
      xx <- xy[,-ncolxy,drop=FALSE]
      KgLow <- apply(xx,2,min)
      KgUp <- apply(xx,2,max)
      maxrange <- KgUp-KgLow
      GCVlowerFactor <- 2 # or 20 ? cf comments on Migraine default in Migraine code
      GCVupperFactor <- 5
      lower <- maxrange/(GCVlowerFactor*nuniquerows)
      upper <- maxrange*GCVupperFactor
      if (maxSmoothness-minSmoothness>1e-6) { # fixed Smoothness case
        maxrange <- c(maxrange,maxSmoothness-minSmoothness)
        lower <- c(lower,minSmoothness)
        upper <- c(upper,maxSmoothness)
        fixedSmoothness <- numeric(0) ## not c() which is NULL...
      } else {
        initCovFnParam <- initCovFnParam[seq_len(length(lower))]
        fixedSmoothness <- maxSmoothness
      }
      for (it in seq_len(length(lower))) {
        if (initCovFnParam[it]<lower[it] ## includes default case
            || initCovFnParam[it]>upper[it]) {
          initCovFnParam[it]=lower[it]+maxrange[it]*2/3;
        }
      }
      if ("L-BFGS-B" %in% method) {
        control <- list(parscale=(upper-lower), trace=verbosityFromoptimizer)
        optr <- optim(par=initCovFnParam,fn=GCV_lamVar_covFix_Wrapper,method="L-BFGS-B",
                      lower=lower,upper=upper,control=control,
                      fixedSmoothness=fixedSmoothness,returnFnvalue=TRUE)
        solution <- optr$par
      } else if ("lbfgsb3" %in% method){ ## very slow
        control <- list(trace=verbosityFromoptimizer)
        if ( ! requireNamespace("lbfgsb3",quietly=TRUE) ) {
          stop("Package lbfgsb3 not installed.")
        }
        optr <- lbfgsb3::lbfgsb3(prm=initCovFnParam,fn=GCV_lamVar_covFix_Wrapper,lower=lower,upper=upper,
                        control=control,
                        fixedSmoothness=fixedSmoothness,returnFnvalue=TRUE)
        solution <- optr$prm
      } else if ("bobyqa" %in% method){ ## marginally better than optim ?
        if ( ! requireNamespace("minqa",quietly=TRUE) ) {
          stop("Package minqa not installed.")
        }
        control <- list(rhobeg=min(abs(upper-lower))/20,iprint=verbosityFromoptimizer)
        control$rhoend <- max(1,control$rhobeg)/1e6
        optr <- minqa::bobyqa(par=initCovFnParam,fn=GCV_lamVar_covFix_Wrapper,lower=lower,upper=upper,
                       control=control,
                       fixedSmoothness=fixedSmoothness,returnFnvalue=TRUE)
        solution <- optr$par
      } else if ("NLOPT_LN_BOBYQA" %in% method){ ## may be the fastest
        optr <- nloptr(x0=unlist(initCovFnParam),eval_f=GCV_lamVar_covFix_Wrapper,lb=lower,ub=upper,
                       opts=list(algorithm="NLOPT_LN_BOBYQA",xtol_rel=1.0e-4,maxeval=-1,print_level=verbosityFromoptimizer),
                       fixedSmoothness=fixedSmoothness,returnFnvalue=TRUE)
        solution <- optr$solution
      } else stop("Unknown 'method' in CKrigcoefs()")
      resu <- list(covfnparam=solution,fnEvalCount=getFnEvalCount(),method=method)
    } else { ## no estim of covfnparam but we will estimate lambda
      resu <- list(covfnparam=covfnparamA)
      fixedSmoothness <- numeric(0) ## local value for GCV_lamVar_covFix_Wrapper() call:
    }
    if (is.na(lambdaA)) resu$lambda <- GCV_lamVar_covFix_Wrapper(resu$covfnparam,fixedSmoothness=fixedSmoothness,returnFnvalue=FALSE)
    deleteCSmooth() ## this code conditional on optimisation
  }
  return(resu) # is a list with elements depending on call:
  # if optimise covfnparam -> covfnparam, lambda, fnEvalCount, method
  # else if optimise lambda -> input covfnparamA, lambda
  # else c("c","d","CKrigidx")
} ## end def CKrigcoefs

# R-style implementation of default values
R_GCV_lamVar_covFix <- function(a,fixedSmoothness=numeric(0L),returnFnvalue=TRUE) {
  ## call to R-to-C wrapper with explicit values
  GCV_lamVar_covFix_Wrapper(a=a,fixedSmoothness=fixedSmoothness,returnFnvalue=returnFnvalue)
}
