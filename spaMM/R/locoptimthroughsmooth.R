sampleNearby <- function(focalPts,n=NULL,stepsizes,
                         margin ## to avoid too close pairs of points
                         ) {
  d <- ncol(focalPts)
  nr <- nrow(focalPts)
  if (n<nr) {
    subfocal <- focalPts[sample(seq_len(nr),n),,drop=FALSE]
  } else subfocal <- focalPts
  rmax <- stepsizes/2
  rmin <- margin*rmax
  randpts <- apply(subfocal,1,function(v) {
    rdxy <- sample(c(-1,1),d,replace=TRUE)* runif(d,min=rmin,max=rmax) 
    v + rdxy 
  })
  if (d==1L) {
    randpts <- as.matrix(randpts) 
  } else randpts <- t(randpts)
  colnames(randpts) <- colnames(focalPts) ## (in 1D at least) rbind checks that the names match...
  randpts
}

## uses gridSteps or n depending on sampling
## default sampling depends on length(lower)
sampleGridFromLowUp <- function(LowUp,n=NULL,gridSteps=NULL,sampling=NULL) {
  ## grid from lower, upper
  lower <- LowUp$lower
  upper <- LowUp$upper
  d <- length(lower)
  if (is.null(sampling)) {
    if (d<4) {sampling="grid"} else {sampling="rvolTriangulation"}
  }
  byvar <- t(rbind(unlist(lower),unlist(upper))) 
  byvar <- 0.999 * byvar + 0.001 *rowMeans(byvar)
  grillelist <- list()
  if (is.null(gridSteps)) {
    gridSteps <- c(10,6,4,3)[min(d,4)] ## => 10  36  64  81 243 729 points 
  }
  for(name in rownames(byvar)) {grillelist[[name]] <- seq(byvar[name,1],byvar[name,2],length.out=gridSteps)}
  pargrid <- expand.grid(grillelist)
  if (sampling=="rvolTriangulation") { ## while randomly sample the simplices defined by the regular grid
    vT <- volTriangulation(pargrid) ## note no convhulln call -> many internal simplices
    if (is.null(n)) n <- gridSteps^min(d,6) ## default maxi 729 for high d
    pargrid <- rbind(pargrid,rvolTriangulation(n,vT)) ## regular + random
  } else { ## more structured sampling: will sample the grid 
    attr(pargrid,"regularGrid") <- seq_len(nrow(pargrid))
    ## redefines a grid
    insides <- lapply(grillelist, function(v) {(v[2]-v[1])/2+v[-c(gridSteps)]}) ## 19 5 3 2 2 2...
    stepsizes <- unlist(lapply(insides, function(v) {v[2]-v[1]})) 
    insides <- expand.grid(insides) ## 9 25 27 16 32 64 ...
    randgrid <- sampleNearby(insides,n=nrow(insides),stepsizes=stepsizes,margin=1e-4) ## not really nearby given the large stepsizes
    pargrid <- rbind(pargrid,randgrid) ## regular + random = 19 61 91 97 275... (20 replicates typically added)
  }
  ##
  pargrid
}

#pargrid <- gridFromLowUp(LowUp)

# doc provisoire:
# pargrid : a matrix of fittedPars
# anyHLCor.args = arguments for evaluation of likelihood by HLCor
# prevPtls : previous points (fittedPars, respName, lambda)
# control.smooth : smoothing parameters (typically $rho )and/or number of duplicates ($nrepl)
#
# This function estimates likelihood in the input points + a second estimate in nrepl of these points 
# The computes a smoothed likelihood surface (default $nu=4) by ordinary kriging
# Then finds the maximum of this smoothed likelihood surface
#
# retruns a list with elements
# $par: cf optimize
# $value: cf optimize
# $predictions: predicted responses values at smoothed point values; attributes include an estimate of the variance of [likelihood estimation by HLCor]
# $Krigobj: ... 
# $forSmooth: the input data of the smoothing computation, i.e. the estimates of likelihood by HLCor
#
#

NAcleaning <- function(forSmooth,inPairs) {
  isNAseInt <- is.na(forSmooth$seInt) ## after reordering !
  misleadingPairs <- inPairs & (isNAseInt[-length(isNAseInt)] | isNAseInt[-1]) 
  misleadingPairs <- c(FALSE,misleadingPairs) | c(misleadingPairs,FALSE) ## pairs with at least one NA 
  NAcleaned <- forSmooth[ ! misleadingPairs,,drop=FALSE]
  # NA may mean we have no info on seInt in which case we can fit phi using seInt
  #  and then we should ideally check whether resid formula uses seInt to decide whether NA singletons should be removed 
  #   (or to know whether they will automatically be removed by HLfit ?)
  # but NA may have a more general meaning... let's be cautious
  NAcleaned <- NAcleaned[ ! is.na(NAcleaned$seInt),,drop=FALSE] ## removes remaining NA's not from pairs
  return(NAcleaned)
}

optimthroughSmooth <- function(pargrid,anyHLCor.args,prevPtls=NULL,control.smooth=list(),verbose=interactive()) {
  ranges <- apply(pargrid,2,range)
  LowUp <- apply(ranges,1,as.list)
  names(LowUp) <- c("lower","upper")
  lower <- LowUp$lower
  upper <- LowUp$upper
  ## 
  processedHL1 <- getProcessed(anyHLCor.args$processed,"HL[1]",from=1L) ## there's also HLmethod in processed<[[]]>$callargs
  logLobj <- anyHLCor.args$`HLCor.obj.value`
  prevmsglength <- 0
  ## eval nrepl
  nrepl <- control.smooth$nrepl
  if (processedHL1 == "SEM") {
    if (is.null(nrepl)) {
      nrepl <- Inf ## default is to duplicate all simuls !
    }
    ## don't confuse 'lower' length (corresponding to the rho,nu of the spatial process): may be 1 or 2 depending eg on user's ranFix$nu  
    ## an length of control.smooth$ranFix (for likelihood surface smoothing) ## which should be 2 (for fixed nu) if no smoothing params are to be reestimated
    ## but length(control.smooth$ranFix) can meaningfully be compared to length (lower). Note that $ranFix will typically have rho, nu while lower has transformed params.
    if (nrepl == 0 && length(control.smooth$ranFix$rho)!=length(lower)) {
      message("(!) From optimthroughSmooth: suspect nrepl==0 for stochastic simulation\n  with likelihood surface correlation parameters to be estimated.")
    }
  } else {
    if (is.null(nrepl)) nrepl <- 0
    if (nrepl>0) message("(!) From optimthroughSmooth: suspect nrepl>0 for deterministic simulation.")
  }  
  #
  subpargrid <- pargrid[sample(nrow(pargrid),min(nrepl,nrow(pargrid))),,drop=FALSE]
  pargrid <- as.matrix(rbind(pargrid, subpargrid))
  grid.obj <- matrix(NA,nrow=NROW(pargrid),ncol=4L)
  colnames(grid.obj) <- c("logLobj","lambda","seInt","pmvnorm")
  if (verbose) cat("\n") 
  skeleton <- anyHLCor.args$skeleton
  for (ii in seq_len(NROW(pargrid))) {
    anyHLCor.args$ranPars[names(skeleton)] <- relist(pargrid[ii,],skeleton)
    ## as in HLCor.obj:
    attr(anyHLCor.args$ranPars,"RHOMAX") <- attr(skeleton,"RHOMAX")
    attr(anyHLCor.args$ranPars,"NUMAX") <- attr(skeleton,"NUMAX")
    types <- attr(skeleton,"type")
    attr(anyHLCor.args$ranPars,"type")[names(types)] <- types
    hlcor <- do.call("HLCor",anyHLCor.args) ## this reconstructs a list of the form of initvec, then adds it to other anyHLCor$ranPars information
    vec <- c("logLobj"=hlcor$APHLs[[logLobj]],lambda=hlcor$lambda,
      seInt =attr(hlcor$APHLs[[logLobj]],"seInt"), ## seInt attr may be NULL then no seInt element
      pmvnorm = grepl("pmvnorm",attr(hlcor$APHLs[[logLobj]],"method"))
      ) ## seInt attr may be NULL then no seInt element
    grid.obj[ii,seq_len(length(vec))] <- vec
    #     if (is.infinite(grid.obj[ii,"logLobj"])) { ## if underflow in estim of L (eg, pmvnorm)
    #       # .Machine$double.xmin gives "the smallest non-zero normalized floating-point number" but smaller "denormalized" numbers are possible
    #       # see Double-precision values in ?double 
    #       grid.obj[ii,"logLobj"] <- - 745 ## log(4.940656e-324); a tooLow mechanism as in Infusion may be a better approach
    #       grid.obj[ii,"seInt"] <- NA ## tag for removal from any fit using variable seInt (and a bit more analyses: see below)
    #       ## Fr->FR en fait ce qu'il faut c'est un mecanisme remove.pairswithNas, cf un ecran plus bas ?
    #     }
    if (verbose) {
      msg <- paste(ii,"simulations run out of", NROW(pargrid)," ")
      prevmsglength <- overcat(msg, prevmsglength)
    }
  }
  forSmooth <- data.frame(cbind(pargrid,grid.obj))
  if ( ! is.null(prevPtls)) forSmooth <- rbind(prevPtls,forSmooth) ## older points first!
  RHS <- paste(names(lower),collapse="+")
  form <- as.formula(paste("logLobj ~ 1+Matern(1|",RHS,")"))
  ## ****** predict a (profile) likelihood surface for the correlation parameters ******
  ## phi here represent varSEMandMCint the variance of MC integration and that of the SEM algo
  ##  phi = a [SEM]+ b lambda [MCint]: 
  ## i e log(phi) = log(a+ b lambda), not a +b log(lambda) => code pas ideal
  ## fortunately lambda may be roughly indep of the corrpars
  ## Further the IRLS in glm.fit is fussy... 
  ## et le fit des dev resids peut etre surprenant, cf plot(var explicative,log(resp)) dans dispGammaGLM -> glm(resp ~...)
  ## hence we dont use this residual.formula:
  ###Krigobj <- corrHLfit(form,data=forSmooth,resid.formula= ~log(lambda),init.corrHLfit=list(rho=rep(1,length(lower))),ranFix=list(nu=4))
  ranFix <- list(nu=4)
  if (processedHL1 != "SEM") ranFix$phi <- 1e-06 ## interpolation
  ranFix[names(control.smooth$ranFix)] <- control.smooth$ranFix
  initSmooth <- control.smooth$initSmooth ## corr pars du smoothing
  if (is.null(initSmooth)) initSmooth <- rep(1,length(lower))
  init.corrHLfit <- list(rho=initSmooth) ## important as it gives the length of rho to corrHLfit
  init.corrHLfit[names(ranFix)] <- NULL
  # the data for the phi GLM contain controlled names (trRho trNu logLobj lambda seInt pmvnorm)
  ## defaults:
  resid.formula <- control.smooth[["resid.model"]]$formula
  ## we look whether there is suspect variation in seInt
  if(is.null(resid.formula)) {
    ## Do not use seInt when there are enough pmvnorm, presumably at the top of the lik surf (even though seInt is available)
    if ( is.null(forSmooth$seInt) ) { ## non-standard case
      resid.formula <- ~ 1 
    } else if (length(validSEs <- which(forSmooth[,"seInt"]>0L))>20L) { ## removes NaN's...
      seInt <- forSmooth[validSEs,"seInt"]
      seInt975 <- mean(seInt) + sqrt(var(seInt))* 1.96
      if (forSmooth[which.max(forSmooth[,"logLobj"]),"seInt"] > seInt975) {
        resid.formula <- ~1+offset(seInt^2) 
      } else resid.formula <- ~ 1
    } else resid.formula <- ~ 1 ## many NaN's [, or (oddly), many 0's from pmvnorm, if GHK is not called in that case]
  }
  resid.family <- control.smooth[["resid.model"]]$family ## may be NULL; not that of the user-level model
  if(is.null(resid.family)) {
    if (deparse(resid.formula[[2]]) == "1") {
      resid.family <- spaMM_Gamma(log)
    } else resid.family <- spaMM_Gamma(identity)  
  }   ## this makes ~1+offset(seInt^2) + identity the default for GHK, but the hyper-default is ~1 with pmvnorm
  resid.model <- list(formula=resid.formula, family=resid.family)
  #
  ## logic: (if some trouble) {
  ##  use subset of data to estimate corrPars , forCorrEst is not NULL
  ## } ## else use all data, single corrHLfit sufficient, forCorrEst is NULL
  if ( ! is.null(forSmooth$seInt) && any(is.na(forSmooth$seInt))) { ## is.na(NaN) is TRUE
    forSmooth <- forSmooth[do.call(order,forSmooth),]
    inPairs <- apply(diff(as.matrix(forSmooth[,1:2,drop=FALSE]))==0,1,all)
    ## cf Infusion:::remove.pairswithNas for detailed explanations of the test
    forCorrEst <- NAcleaning(forSmooth,inPairs) ## NULL if no NAs/NaN to clean
  } else  forCorrEst <- NULL
  if (! is.null(forCorrEst)) {
    message("NA/NaN in SEs of likelihood estimates")
    Krigobj <- corrHLfit(form,data=forCorrEst,resid.model= resid.model,
                         init.corrHLfit=init.corrHLfit,ranFix=ranFix)
    ## now we can use "consistent NAs" in the following fit IF we can provide phi and lambda estimates...
    if (deparse(resid.formula[[length(resid.formula)]]) == "1") {
      inconsistentNAs <- inPairs & (diff(is.na(forSmooth$seInt)) != 0L)
      misleadingPairs <- c(FALSE,inconsistentNAs) | c(inconsistentNAs,FALSE) ## pairs with only one NA 
      forSmooth <- forSmooth[ ! misleadingPairs,,drop=FALSE]
      ranFix$lambda <- Krigobj$lambda
      ranFix$phi <- Krigobj$phi ## *FR->FR ineffective for resid.predictor$formula= ~1 hence phi is refit below*
    } else forSmooth <- forCorrEst ## play safe
    ## in all cases:
    ranFix$rho <- Krigobj$corrPars$rho
  }
  init.corrHLfit[names(ranFix)] <- NULL
  Krigobj <- corrHLfit(form,data=forSmooth,resid.model= resid.model ,
                       init.corrHLfit=init.corrHLfit,ranFix=ranFix)
  ### quick check of variance of logL estimation:
  if (FALSE) {
    essai <- forSmooth[do.call(order,forSmooth[,seq_len(length(lower))]),,drop=FALSE]
    diffs <- apply(essai,2,diff)
    ## next line selects diff lines for identical coordinates and takes from them the diff value for logL 
    logLdiffs <- diffs[apply(diffs[,seq_len(length(lower))]==rep(0,length(lower)),1,all),length(lower)+1]
    var(logLdiffs)/2 ## since diff= sum(eps_1+eps_2)
  }
  ###
  ## ****** optimize in the predicted likelihood surface ******
  predictions <- predict(Krigobj,binding="fitted")
  predictions <- predictions[order(predictions[,attr(predictions,"fittedName")],decreasing=TRUE),]
  initvec <- predictions[1,names(lower)]
  ## redefines lower, upper, for maximization
  ranges <- apply(forSmooth[,names(lower),drop=FALSE],2,range)
  LowUp <- apply(ranges,1,as.list)
  names(LowUp) <- c("lower","upper")
  lower <- LowUp$lower
  upper <- LowUp$upper
  ##
  optr <- optim(initvec,function(v) {predict(Krigobj,v)[,1]},method="L-BFGS-B",
                control=list(fnscale=-1),lower=unlist(lower),upper=unlist(upper))
  names(optr$par) <- names(lower)
  attr(predictions,"fittedPars") <- names(lower) 
  attr(predictions,"MSy") <- Krigobj$phi ## FR->FR ecrire un extracteur pour phi... 
  optr$predictions <- predictions 
  optr$Krigobj <- Krigobj ## 
  optr$forSmooth <- forSmooth
  return(optr)
} ## end def optimthroughSmooth
