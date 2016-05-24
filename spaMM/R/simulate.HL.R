## derived from the devel version of lme4:
##' tests if argument implies a model without random effects or not
##' TRUE if argument is NA or ~0, 
##' FALSE if argument is NULL or a non-trivial formula
noReForm <- function(re.form) {
  (!is.null(re.form) && !is(re.form,"formula") && is.na(re.form)) ||
    (inherits(re.form,"formula") && length(re.form)==2 && identical(re.form[[2]],0))
}

newetaFix <- function(object, newMeanFrames) {
  ## newdata -> offset must be recomputed. 
  off <- model.offset( newMeanFrames$mf) ### look for offset from (ori)Formula 
  if ( is.null(off) ) { ## ## no offset (ori)Formula term. 
    ## then we check that no non zero $offset was used. This would make prediction generally incorrect
    if (! is.null(off <- attr(object$predictor,"offsetObj")$offsetArg)) { ## that means there was a non trivial offset argument in the original Predictor(formula...)  
      message("Prediction in new design points from a fit with formula=Predictor(... non-NULL offset ...) is suspect.")
    } ## but we still proceed with this dubious offset
  }    
  ## dans l'état actuel $fixef et complet,incluant les etaFix$beta: pas besoin de les séparer
  if (ncol(newMeanFrames$X)>0) {
    etaFix <-  drop(newMeanFrames$X %*% object$fixef)
  } else {
    etaFix <- 0 
  }   
  if ( ! is.null(off)) etaFix <- etaFix + off   
  return(etaFix)
}



# simulate.HLfit(fullm[[2]],newdata=fullm[[1]]$data,size=fullm[[1]]$data$total) for multinomial avec binomial nichées de dimension différentes
# FR->FR misses the computation of randoem effects for new spatial positions: cf comments in the code below
simulate.HLfit <- function(object, nsim = 1, seed = NULL, newdata=NULL, sizes=object$weights,...) { ## object must have class HLfit; corr pars are not used, but the ZAL matrix is.
  ## RNG stuff copied from simulate.lm
  if (!exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE))
    runif(1)
  if (is.null(seed))
    RNGstate <- get(".Random.seed", envir = .GlobalEnv)
  else { ## this makes changes to RNG local where 'seed' is used:
    R.seed <- get(".Random.seed", envir = .GlobalEnv)
    set.seed(seed)
    RNGstate <- structure(seed, kind = as.list(RNGkind()))
    on.exit(assign(".Random.seed", R.seed, envir = .GlobalEnv))
  }
  if (inherits(object,"HLfitlist")) { ## testing for list is not valid since an HLfit object is always a list
    message("simulate does not yet work on list of fits as returned by multinomial fit:")
    message(" run simulate on each of the individual fit in the list")
    stop() ## FR->FR also some basic changes in fixedLRT but more would be needed 
  }  
  if (is.null(newdata)) {
    # rebuild linear predictor in three steps eta= X . beta + off + ZAL . RANDOM v
    # hence do not use object$eta which contains PREDICTED V
    nobs <- length(object$y)
    if (ncol(object$`X.pv`)>0) {eta <- object$`X.pv` %*% object$fixef} else {eta <- rep(0,nobs)}
    ##
    eta <- eta + attr(object$predictor,"offsetObj")$total ## a PROCESSED predictor or resid.predictor always has a non-NULL offset term  
  } else {
    nobs <- nrow(newdata)
    ## [-2] so that HLframes does not try to find the response variables  
    allFrames <- HLframes(formula=attr(object$predictor,"oriFormula")[-2],data=newdata) ## may need to reconstruct offset using formula term
    eta <- newetaFix(object,allFrames) ## X . beta + off
  }
  ##
  if (any(object$models[["lambda"]] != "")) { ## i.e. not a GLM
    if (is.null(newdata)) {
      ZAL <- attr(object$predictor,"ZALMatrix")
      cum_n_u_h <- attr(object$lambda,"cum_n_u_h")
      vec_n_u_h <- attr(cum_n_u_h,"vec_n_u_h")
    } else {
      FL <- spMMFactorList(object$predictor, allFrames$mf, 0L, drop=TRUE) 
      ##### simulation given the observed response ! :
      #       spatial.terms <- findSpatial(locform)
      #       spatial.model <- spatial.terms[[1]] 
      #       if( is.null(spatial.model)) {
      #         uuCnewold <- NULL
      #       } else {
      #         v_h_coeffs <- predictionCoeffs(object) ## changes the coefficients in the right u_range
      #         blob <- calcNewCorrs(object=object,locdata=newdata,spatial.model=spatial.model)
      #         uuCnewold <- blob$uuCnewold
      #       }
      #       ZALlist <- computeZAXlist(XMatrix=uuCnewold,ZAlist=FL$Design)
      #       (unfinished:) il faut rajouter la conditional variance comme dans le SEM => simuler comme dans le SEM ?
      ##### independent simulation ! : 
      # en fait il faut recycler du code de HLCor...
      ##### the following code with NULL LMatrix ignores spatial effects with newdata:
      ZALlist <- computeZAXlist(XMatrix=NULL,ZAlist=FL$Design)
      nrand <- length(ZALlist)
      vec_n_u_h <- unlist(lapply(ZALlist,ncol)) ## nb cols each design matrix = nb realizations each ranef
      cum_n_u_h <- cumsum(c(0,vec_n_u_h))
      ZALlist <- lapply(seq_len(length(ZALlist)),as.matrix)
      ZAL <- do.call(cbind,ZALlist)
    }
    lcrandfamfam <- attr(object$rand.families,"lcrandfamfam") ## unlist(lapply(object$rand.families,function(rf) {tolower(rf$family)})) 
    fittedLambda <- object$lambda.object$lambda_est
    newV <- lapply(seq(length(vec_n_u_h)), function(it) {
      nr <- vec_n_u_h[it]
      u.range <- (cum_n_u_h[it]+1L):(cum_n_u_h[it+1L])
      loclambda <- fittedLambda[u.range]
      newU <- replicate(nsim,{
        switch(lcrandfamfam[it], ## remainder of code should be OK for rand.families
               gaussian = rnorm(nr,sd=sqrt(loclambda)),
               gamma = rgamma(nr,shape=1/loclambda,scale=loclambda),
               beta = rbeta(nr,1/(2*loclambda),1/(2*loclambda)),
               "inverse.gamma" = 1/rgamma(nr,shape=1+1/loclambda,scale=loclambda), ## yields inverse gamma (1+1/object$lambda,1/object$lambda)
               stop("(!) random sample from given rand.family not yet implemented")
        )},simplify=TRUE) ## should have nsim columns
      object$rand.families[[it]]$linkfun(newU) 
    }) ## one multi-rand.family simulation
    newV <- do.call(rbind,newV) ## each column a simulation
    if (nsim==1L) {
      eta <- eta + ZAL %id*% newV 
    } else eta <-  matrix(rep(eta,nsim),ncol=nsim) + ZAL %id*% newV ## nobs rows, nsim col
  }
  mu <- object$family$linkinv(eta) ## ! freqs for binomial, counts for poisson
  ## 
  phiW <- object$phi/object$prior.weights ## cf syntax and meaning in Gamma()$simulate / spaMM_Gamma$simfun
  famfam <- tolower(object$family$family)
#  if (famfam=="binomial" && is.null(size)) size <- object$weights ## original sample size by default
  respv <- function(mu) {switch(famfam,
                                gaussian = rnorm(nobs,mean=mu,sd=sqrt(phiW)),
                                poisson = rpois(nobs,mu),
                                binomial = rbinom(nobs,size=sizes,prob=mu),
                                gamma = rgamma(nobs,shape= mu^2 / phiW, scale=phiW/mu), ## ie shape increase with prior weights, consistent with Gamma()$simulate / spaMM_Gamma()$simfun
                                stop("(!) random sample from given family not yet implemented")
  )} ## vector
  if (nsim>1) {
    resu <- apply(mu,2,respv) ## matrix
  } else {
    resu <- respv(mu)
  }
  return(resu)    
}


simulate.HLfitlist <- function(object,nsim=1,seed=NULL,newdata=object[[1]]$data,sizes=object[[1]]$weights,...) {
  ## RNG stuff copied from simulate.lm
  if (!exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE))
    runif(1)
  if (is.null(seed))
    RNGstate <- get(".Random.seed", envir = .GlobalEnv)
  else { ## this makes changes to RNG local where 'seed' is used:
    R.seed <- get(".Random.seed", envir = .GlobalEnv)
    set.seed(seed)
    RNGstate <- structure(seed, kind = as.list(RNGkind()))
    on.exit(assign(".Random.seed", R.seed, envir = .GlobalEnv))
  }
  replicate(nsim, {
    allrownames <- unique(unlist(lapply(object,function(hl){rownames(hl$data)})))
    resu <- matrix(0,nrow=length(allrownames),ncol=length(object)) ## two cols if 3 types
    cumul <- 0
    if (length(sizes) != nrow(newdata)) {
      mess <- pastefrom("length(sizes) != nrow(newdata).",prefix="(!) From ")
      stop(mess)
    }
    for (it in seq(ncol(resu))) {
      ## it = 1 se ramène à simulate(object[[1]])
      resu[,it] <- simulate(object[[it]],newdata=newdata,sizes=sizes - cumul)
      cumul <- rowSums(resu)  
    }
    resu <- cbind(resu,sizes - cumul) ## now 3 cols if 3 types
    rownames(resu) <- allrownames
    colnames(resu) <- attr(object,"sortedTypes")
    as.data.frame(resu)
  },simplify=FALSE)
}
  

## there is update.HL for new fits of the same X and same Y...

## there is ?? for new X and new Y... GCV vs Krig....






