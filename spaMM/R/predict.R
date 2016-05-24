`%id*%` <- function(A,b) {
  if (is.identity(A)) return(b)
  return(A %*% b)
}


`%id*id%` <- function(A,B) {
  if (is.identity(A)) return(B)
  if (is.identity(B)) return(A)
  return(A %*% B)
}

calc_beta_w_cov <- function(object) {
  beta_w_cov <- attr(object$beta_cov,"beta_v_cov")
  invL <- calc_invL(object) ## correlation matrix of ranefs is solve((t(invL)%*%(invL)))
  # invL must be solve(attr(object$predictor,"LMatrix"))
  if ( ! is.null(invL)) {
    pforpv <- ncol(object$beta_cov)
    v.range <- pforpv+seq(ncol(invL))
    beta_w_cov[v.range,] <- t(invL) %*% beta_w_cov[v.range,]
    beta_w_cov[,v.range] <- beta_w_cov[,v.range] %*% invL
  }
  return(beta_w_cov)
}

calcPredVar <- function(Coldnew,X.pv,newZAC=NULL,newZA,beta_w_cov,
                        CnewnewList=NULL,
                        invColdoldList=NULL, ## may be null is no newdata.
                        lambda=NULL, ## may be null if no newdata.
                        logdispObject=NULL, ## should remain NULL is disp not requested
                        covMatrix=FALSE,blockSize=100L) {
  nrX <-  nrow(X.pv)
  if (is.matrix(Coldnew) || inherits(Coldnew,"Matrix") || inherits(Coldnew,"ff")) {
    Coldnew <- list(dummyid=Coldnew)
  } ## if (!is.list(Coldnew)) not valid bc ff objects inherit from list()
  if (is.matrix(newZA) || inherits(newZA,"Matrix") || inherits(newZA,"ff")) {
    newZA <- list(dummyid=newZA)
  }
  if ( ( ! covMatrix ) && nrX > blockSize) {
    slices <- unique(c(seq(0L,nrX,blockSize),nrX))
    sliceVar <- function(it) {
      slice <- (slices[it]+1L):slices[it+1L]
      ## here  the problem is that newZA should map the new response levels 
      ## to the 'new' levels of random effects  
      nrand <- length(newZA) ## or of any other of the lists of matrices
      requiredLevelsList <- lapply(seq_len(nrand), function(rd) {
        which(colSums(newZA[[rd]][slice,,drop=FALSE])>0L) })
      locnewZA <- lapply(seq_len(nrand), function(rd) { 
        newZA[[rd]][slice,requiredLevelsList[[rd]],drop=FALSE] })
      # predVar in observed points uses C rather than L hence we need to compute ZA.C in all cases ('newZAC")
      #  andthen we need newZA and Coldnew
      locColdnew <- lapply(seq_len(nrand), function(rd) { 
        locC_on <- Coldnew[[rd]][,requiredLevelsList[[rd]],drop=FALSE]
        attr(locC_on,"isEachNewInOld") <- attr(Coldnew[[rd]],"isEachNewInOld")[requiredLevelsList[[rd]]] ## for non-spatial effects; (qualifies sub cols of sub Cnewold)
        return(locC_on)
      })
      locnewZAClist <- lapply(seq_len(nrand), function(rd) { 
        terme <- locnewZA[[rd]] %id*id% t(locColdnew[[rd]])[] 
        as.matrix(terme) ## loses names but they are not useful here 
      })
      if (nrand>1L) {newZAC <- do.call(cbind,locnewZAClist)} else {newZAC <- locnewZAClist[[1]]}
      ## the only one that corresponds to real newdata in predict()
      if (!is.null(CnewnewList)) {
        locCnewnewList <- lapply(seq_len(nrand), function(rd) { 
          CnewnewList[[rd]][requiredLevelsList[[rd]],requiredLevelsList[[rd]],drop=FALSE] }) 
      } else locCnewnewList <- NULL
      calcPredVar(Coldnew=locColdnew, ## needed only if newdata in predict() call
                  X.pv=X.pv[slice,,drop=FALSE],## problem is that this creates the apparence of new data end more calculations
                  newZAC=newZAC,
                  newZA=locnewZA, ## either newZA or newZAC needed even if no newdata in predict() call
                  beta_w_cov=beta_w_cov,CnewnewList=locCnewnewList,
                  invColdoldList=invColdoldList, ## needed only if newdata in predict() call
                  lambda=lambda,
                  logdispObject=logdispObject,
                  covMatrix=covMatrix,blockSize=blockSize)
    }
    unlist(sapply(seq_len(length(slices)-1L),sliceVar))
  }
  ############################################################
  if (is.null(newZAC)) {
    nrand <- length(newZA) ## or of any other of the lists of matrices
    newZAClist <- lapply(seq_len(nrand), function(it) {
      terme <- newZA[[it]] %id*id% t(Coldnew[[it]])[] ## %id*id% further keeps col names of C if tZA==I
      as.matrix(terme) ## loses names but they are not useful here 
    })
    if (nrand>1L) {newZAC <- do.call(cbind,newZAClist)} else {newZAC <- newZAClist[[1]]}
  }
  newAugX <- CBIND(X.pv,newZAC) ## mais en fait pas un AugX since it uses C (in C.w) rather than L (in L.v)
  ## First component of predVar
  # variance of expectation of Xbeta+Zb due to var of (hat(beta),hat(v)) using E[b] as function of hat(v)
  calcZWZt_mat_or_diag <- function(Z,W,returnMat) { ## fixefVar or fixefVar + a bit of ranefVar
    if (returnMat) {
      return(Z[] %id*id% W[] %id*id% t(Z)[])
    } else {
      premul <- Z[] %id*id% W[]
      return(rowSums(suppressMessages(premul * Z[]))) ## suppress message("method with signature...") [found by debug(message)]
    }
  }
  predVar <- calcZWZt_mat_or_diag(newAugX,beta_w_cov,covMatrix)
  ## Second component of predVar:
  # Evar: expect over distrib of (hat(beta),hat(v)) of variance of Xbeta+Zb given (hat(beta),hat(v))
  if (! is.null(CnewnewList) ) {
    nrand <- length(newZA) ## or of any other of the lists of matrices
    Evarlist <- lapply(seq_len(nrand), function(it) {
      isEachNewInOld <- attr(Coldnew[[it]],"isEachNewInOld")
      if ( ! is.null(isEachNewInOld)) { ## non spatial effect: a vector of booleans indicating whether new is in old (qualifies cols of Cnewold)
        Evar <- Diagonal(x=lambda[it]*as.numeric(! isEachNewInOld))
      } else { ## spatial effect
        Cno_InvCoo_Con <- t(Coldnew[[it]])[] %id*id% (invColdoldList[[it]])[] %id*id% (Coldnew[[it]])[]
        Evar <- lambda[it] * (CnewnewList[[it]] - Cno_InvCoo_Con)
      } 
      terme <- calcZWZt_mat_or_diag(newZA[[it]],Evar,covMatrix)
      if (covMatrix) {
        return(as.matrix(terme))  ## loses names but they are not useful here
      } else return(terme)
    })
    if (nrand>1L) {Evar <- Reduce("+",Evarlist)} else {Evar <- Evarlist[[1]]}
    predVar <- predVar + Evar
  } 
  if (! is.null(logdispObject)) { ## ie if disp was requested
    # scan problems
    problems <- logdispObject$problems
    if (length(problems)>0L) {
      warning("Some problems were encountered, affecting computation of prediction variance component\n for uncertainty in dispersion parameters:")
      lapply(problems,warning)
    }
    #
    newZACw <- newZAC %*% logdispObject$dwdlogdisp ## typically (nnew * n_u_h) %*% (n_u_h * 2) = nnew * 2 hence small 
    if (covMatrix) {
      disp_effect_on_newZACw <- newZACw %*% logdispObject$logdisp_cov %*% t(newZACw)  
    } else {
      premul <- newZACw %*% logdispObject$logdisp_cov
      disp_effect_on_newZACw <- rowSums(premul * newZACw)
    }
    predVar <- predVar + disp_effect_on_newZACw
  }
  return(predVar) ## may be a Matrix
}



calcResidVar <- function(object,newdata=NULL) {
  phi.object <- object$phi.object
  if (is.null(phi_outer <- phi.object$phi_outer)) { ## valid whether newdata are NULL or not:
    residVar <- predict(phi.object$glm_phi, newdata=newdata, type="response")
  } else { ## phi, but not glm_phi
    if (length(phi_outer)==1L) {
      if (is.null(newdata)) {
        residVar <- rep(phi_outer,nrow(object$X.pv)) ## assumes (length(phi_outer)==1L)           
      } else residVar <- rep(phi_outer,nrow(newdata))
    } else stop("Unable to compute 'residVar' given length(phi_outer)!=1L.") ## and no glm_phi
    ## FR->FR we could improve this if we get a glm_phi when phi was estimed by outer iterations
  }
  residVar
}  

calcNewCorrs <- function(object,locdata,which,
                         spatial.model) {
  resu <- list(uuCnewold=NULL,uuCnewnew=NULL)
  if (any(unlist(which))) {
    # (1) code  common to different ranef models 
    olduniqueGeo <- attr(object,"info.uniqueGeo")
    geonames <- colnames(olduniqueGeo)
    newuniqueGeo <- calcUniqueGeo(data=locdata[,geonames,drop=FALSE]) ## It is essential that it keeps the same order as spMMfactorlist -> ULI -> unique. 
    ### rho only used to compute scaled distances
    rho <- getPar(object$ranFix,"rho")
    #if( !is.null(rho_mapping <- attr(object,"msd.arglist")$rho.mapping)) rho <- rho[rho_mapping] 
    if ( ! is.null(rho_mapping <- attr(object,"msd.arglist")$rho.mapping)
        && length(rho)>1L ) rho <- fullrho(rho=rho,coordinates=geonames,rho_mapping=rho_mapping)
    ## rows from newuniqueGeo, cols from olduniqueGeo:
    msd.arglist <- list(uniqueGeo=newuniqueGeo,uniqueGeo2=olduniqueGeo,
                        rho=rho,return_matrix=TRUE)
    if ( ! is.null(dist.method <- object$control.dist$dist.method)) msd.arglist$dist.method <- dist.method
    if (which$no) resu$uuCnewold <- do.call(make_scaled_dist,msd.arglist) ## ultimately allows products with Matrix ## '*cross*dist' has few methods, not even as.matrix
    if (which$nn)  {
      msd.arglist$uniqueGeo2 <- NULL
      if (nrow(msd.arglist$uniqueGeo)==1L) {
        resu$uuCnewnew <- matrix(0)
      } else resu$uuCnewnew <- do.call(make_scaled_dist,msd.arglist) 
    }
    ## distance matrix and then call to correl fn:
    # (2) code specific to each ranef model
    if ( ! is.null(spatial.model)) {
      corr.model <- as.character(spatial.model[[1]])
      if (corr.model=="AR1") {
        args <-object$ranFix[which(names(object$ranFix) %in% c("ARphi"))]
        if (which$no) resu$uuCnewold <- args$ARphi^resu$uuCnewold  
        if (which$nn) resu$uuCnewnew <- args$ARphi^resu$uuCnewnew  
      } else {
        args <-object$ranFix[which(names(object$ranFix) %in% c("nu","Nugget"))] ## so that rho=1 in MaternCorr
        if (which$no) resu$uuCnewold <- do.call(MaternCorr,args=c(args,list(d=resu$uuCnewold)))  
        if (which$nn) resu$uuCnewnew <- do.call(MaternCorr,args=c(args,list(d=resu$uuCnewnew)))  
      }
    }
  }
  return(resu)
}


makenewname <- function(base,varnames) { ## post CRAN 1.4.1
  varnames <- varnames[which(substring(varnames,1,nchar(base))==base)] 
  allremainders <- substring(varnames,nchar(base)+1) 
  allnumericremainders <- as.numeric(allremainders[which( ! is.na(as.numeric(allremainders )))  ]) ## as.numeric("...")
  ## 2015/03/04
  if (length(allremainders) == 0L && length(allnumericremainders) == 0L) { ## if base = allremainders => length(allnumericremainders) == 0 not sufficient
    fvname <- base
  } else {
    if (length(allnumericremainders)>0L) {
      num <- max(allnumericremainders)+1L
    } else num <- 1L
    fvname <-paste ( base , num , sep="") 
  }
  fvname
}

## (1) for surface prediction: (developed in InferentialSimulation/InferentialSimulation.R)
## (2) But also for generation of fixed effects in simulation of nested-effects models
predict.HLfit <- function(object, newdata = newX, newX = NULL, re.form = NULL,
                          variances=list(),
                          binding = FALSE, 
                          intervals = NULL,
                          level = 0.95,
                          ...) { ## but not new Y
  if ( ! is.null(variances$ranef)) {
    message("'variances$ranef' is obsolete: I interpret this as 'variances$linPred'.")
    variances$linPred <- variances$ranef
    variances$ranef <- NULL
  }
  if ( ! is.null(variances$sum)) {
    message("'variances$sum' is obsolete: I interpret this as 'variances$respVar'.")
    variances$respVar <- variances$sum
    variances$sum <- NULL
  }
  ## the final components returned as attributes have names ...Var, other terms should be named differently
  checkIntervals <- (substr(x=intervals, nchar(intervals)-2, nchar(intervals))=="Var")
  if (any(!checkIntervals)) warning("Element(s)",intervals[!checkIntervals],"are suspect, not ending in 'Var'.")
  # possible elements in return value: fixefVar, predVar, residVar, respVar
  variances[intervals] <- TRUE 
  # $respVar implies components
  if (is.null(variances$predVar)) variances$predVar <- variances$respVar ## may still be NULL
  if (is.null(variances$residVar)) variances$residVar <- variances$respVar ## may still be NULL
  if (is.null(variances$respVar)) variances$respVar <- FALSE 
  # $predVar implies components
  if (is.null(variances$linPred)) variances$linPred <- variances$predVar ## may still be NULL
  if (is.null(variances$disp)) variances$disp <- variances$predVar ## may still be NULL
  if (is.null(variances$predVar)) variances$predVar <- FALSE 
  # Do not let any component empty
  if (is.null(variances$fixefVar)) variances$fixefVar <- FALSE 
  if (is.null(variances$linPred)) variances$linPred <- FALSE 
  if (is.null(variances$disp)) variances$disp <- FALSE ## uncertaintly on dispersion parameters
  if (is.null(variances$residVar)) variances$residVar <- FALSE ## uncertaintly on dispersion parameters
  if (is.null(variances$cov)) variances$cov <- FALSE
  ##
  locform <- attr(object$predictor,"oriFormula")
  ## possible change of random effect terms
  if (noReForm(re.form)) {
    locform <- nobarsMM(locform)  ## ie. (if re.form implies that there is no random effect)
  } else if (inherits(re.form,"formula")) locform <- re.form
  # checking variables in the data
  if (length(locform)==3L) locform <- locform[-2] ## removes RHS for checking  vars on RHS etc
  allvars <- all.vars(locform) 
  if (variances$residVar) allvars <- unique(c(allvars,all.vars(attr(object$resid.predictor,"oriFormula")))) ## but the 'newdata' may not contain the resid.predictor vars. 
  if (is.vector(newdata)) { ## ## less well controlled case, but useful for maximization
    binding <- binding ## :-) binding must be evaluated before newdata is modified
    newdata <- data.frame(matrix(newdata,nrow=1))
    if (length(allvars)==ncol(newdata)) {
      names(newdata) <- allvars
    } else {
      stop(paste("(!) newdata has incorrect length. It should match the following variables:\n",paste(allvars,collapse=" ")))
    }
  } 
  ## it is important that newdata remains NULL if it was so initially because it is tested below. Hence copy in locdata
  if (is.null(newdata)) {
    locdata <- object$data[,allvars,drop=FALSE]
  } else {
    if( is.matrix(newdata) ) newdata <- as.data.frame(newdata)  
    # so that matrix 'newdata' arguments can be used as in some other predict methods.
    locdata <- newdata ## FR->FR [,allvars,drop=FALSE] ?
    if (any(is.na(locdata))) {
      stop("NA's in required variables from 'newdata'. Prediction not possible.")
    }
  }
  ## preparation for fixed effects
  allFrames <- HLframes(formula=locform,data=locdata,fitobject=object)
  newX.pv <- allFrames$X
  if ( ! is.null(betaFix <- attr(object$predictor,"offsetObj")$betaFix)) { ## suppress betaFix cols so that this is consistent with <object>$X.pv 
    newX.pv <- newX.pv[,colnames(object$`X.pv`),drop=FALSE]
  }
  needNewEta <- ( ( ! is.null(newdata) ) || variances$linPred || ! is.null(re.form))
  if (needNewEta) etaFix <- newetaFix(object,allFrames)  ## new fixed effects (or [part of] new linear predictor if re.form)      
  ## preparation for random effects
  spatial.terms <- findSpatial(locform) ## list of formula terms
  spatial.model <- spatial.terms[[1]] ## one formula term, e.g Matern(1|...)
  if ( ! is.null(spatial.model)) { 
    if (! is.null(newdata) && as.character(spatial.model[[1]]) %in% c("adjacency","ar1")) {
      stop("Prediction in newdata not implemented or not possible in the 'adjacency' model")
    } ## FR->FR would be possible for new non-spatial predictor values in the original locations... il faudrait un test sur les elements de la distance matrix
  } 
  ## matching ranef terms of re.form
  if (noReForm(re.form)) {
    nrand <- 0L
  } else {
    if (inherits(re.form,"formula")) {
      ranefs <- parseBars(locform)
      nrand <- length(ranefs)
      newinold <- sapply(ranefs, function(v) {which(v==attr(object$ZAlist,"ranefs"))})
      subZAlist <- object$ZAlist[newinold] ## and reordered
      sublambda <- object$lambda[newinold]
    } else {
      ranefs <- attr(object$ZAlist,"ranefs")
      nrand <- length(ranefs)
      newinold <- seq(nrand)
      subZAlist <- object$ZAlist
      sublambda <- object$lambda
    }    
    spatialOne <- which(ranefs == spatial.model) ## strictly a spatial one, not other correlated ones      
  }
  #
  if (nrand>0L) {
    FL <- spMMFactorList(locform, allFrames$mf, 0L, drop=TRUE) 
    newZAlist <- FL$Design ## must be ordered as parseBars result...
    attr(newZAlist,"ranefs") <- ranefs ## required pour computeZAXlist to match the ranefs of LMatrix
    if ( is.null(spatial.model)) {
      newZAClist <- computeZAXlist(XMatrix=NULL,ZAlist=newZAlist) 
      uuCnewnew <- NULL 
    } else {
      oldLMatrix <- attr(object$predictor,"LMatrix") ## may be NULL
      which <- list(no=( variances$linPred  || needNewEta || ! is.null(newdata)), 
                    nn= variances$linPred && ! is.null(newdata))
      uunewCorrs <- calcNewCorrs(object=object,locdata=locdata,
                                 which=which,
                                 spatial.model=spatial.model)
      ## matrices, not list of matrices which are constructed later
      uuCnewnew <- uunewCorrs$uuCnewnew
      uuCnewold <- uunewCorrs$uuCnewold
      if ( ! is.null(uuCnewold)) attr(uuCnewold,"ranefs") <- attr(oldLMatrix,"ranefs") ## required pour computeZAXlist to match the ranefs of LMatrix
      newZAClist <- computeZAXlist(XMatrix=uuCnewold,ZAlist=newZAlist) ## ZAL's for ZA's and L's (typically some ZA's are unaffected)
    } 
  } 
  #
  ## (1) computes fv (2) compute predVar
  ##### fv
  if (noReForm(re.form)) {
    fv <- object$family$linkinv(etaFix) 
  } else if ( is.null(newdata) && ! inherits(re.form,"formula")) {
    fv <- object$fv ## same length including replicates
    newZAlist <- subZAlist ## useful if predVar
  } else { ## 
    if ( nrand==0L ) {
      eta <- etaFix
      newZAlist <- NULL
    } else {
      #### precomputation of coeffs
      ## on the gaussian scale, L.v_ori ~ lam C (lam C + phi I)^{-1}y 
      ## new random autocorr term ~ lam c (lam C + phi I)^{-1}y = c C^{-1} L_ori.v_ori = c [t(L_ori)]^{-1} v_ori
      ## [t(L_ori)]^{-1} v_ori can be computed once for all predictions => 'w_h_coeffs'
      w_h_coeffs <- object$get_w_h_coeffs() ## should work for nonspatial models 
      old_cum_n_u_h <- attr(object$lambda,"cum_n_u_h")
      lcrandfamfam <- attr(object$`rand.families`,"lcrandfamfam")[newinold]
      augm_w_h_coeffs <- lapply(seq_len(nrand),function(it) {
        oldu.range <- (old_cum_n_u_h[newinold[it]]+1L):(old_cum_n_u_h[newinold[it]+1L])
        if (it %in%  spatialOne) {     # %in% handles zero-length spatialOne...
          return(w_h_coeffs[oldu.range])          
        } else {
          oldlevels <- colnames(subZAlist[[it]])
          newlevels <- colnames(newZAClist[[it]])
          interlevels <- intersect(oldlevels,newlevels)
          oldv <- w_h_coeffs[oldu.range]
          names(oldv) <- oldlevels
          psi_M <- switch(lcrandfamfam[it], 
                           gaussian = 0,
                           gamma = 1, 
                           beta = 1/2, 
                           "inverse.gamma" = 1
          )
          vpsi_M <- object$rand.families[[newinold[it]]]$linkfun(psi_M) 
          ## since vpsi_M can be non-zero, the expectation of the response can be modified in a re.form model compared to the original
          newv <- rep(vpsi_M,length(newlevels))
          names(newv) <- newlevels
          newv[interlevels] <- oldv[interlevels] 
          return(newv)
        }
      })
      if (nrand>1L) {
        ZACw <- lapply(seq_len(nrand),function(it) {
          drop(newZAClist[[it]][] %*% augm_w_h_coeffs[[it]]) ## sapply preforms a cbind if everything is matrix (not Matrix)
        }) ## not sapply which binds 1*1 matrices into a vector of length nrand  
        ZACw <- do.call(cbind,ZACw)
        ZACw <- rowSums(ZACw)
      } else ZACw <- drop(newZAClist[[1]][] %*% augm_w_h_coeffs[[1]])
      eta <- etaFix + ZACw ## (length(eta)) col vector from coeffs = length(eta) row vector...
    }
    # done with eta
    fv <- object$family$linkinv(eta) ## ! freqs for binomial, counts for poisson
  }
  resu <- as.matrix(fv) ## suitable for objective function of optim() etc ## matrix ! maybe more suitable than data frame as objective function
  if ( ! is.logical(binding) ) { ## expecting a string
    binding <- makenewname(base=binding,varnames=colnames(locdata)) ## 09/11/2014 = posterior to CRAN 1.4.1 
    resu <- data.frame(resu)
    colnames(resu) <- binding
    resu <- cbind(locdata,resu) ## becomes a data frame !
    attr(resu,"fittedName") <- binding
  } else { ## expecting binding= FALSE
    if (ncol(locdata)>0)  attr(resu,"frame") <- locdata 
  }
  ##### (2) predVar
  if(variances$linPred) {
    beta_w_cov <- object$get_beta_w_cov()
    if ( ! is.null(newdata)) {
      invColdoldList <- object$get_invColdoldList()
      ## list for Cnewnew, which enters in  newZA %*% Cnewnew %*% tnewZA, hence should not represent newZA itself 
      if (nrand>0L) newnewClist <- lapply(seq_len(nrand),function(it) {
        if ( it %in% spatialOne) { ## avoids running the next algo which is slow on large matrices
          Cnn <- uuCnewnew ## already computed for point prediction
        } else { 
          Cnn <- Diagonal(ncol(newZAlist[[it]])) ## diag(ncol(newZAlist[[it]]))
        }
        return(Cnn)
      }) ## => completely ignores effects removed in re.form 
    } else {
      invColdolList <- NULL
      newnewClist <- NULL
    }
    if (nrand>0L) {
      ## list for Coldnew, which enters in ZA %id*id% Coldnew %id*id% tnewZA 
      oldnewClist <- lapply(seq_len(nrand),function(it) {
        if (it %in% spatialOne) { ## avoids running the next algo which is slow on large matrices
          oldnewC <- t(uuCnewold)
        } else {
          oldlevels <- colnames(subZAlist[[it]])
          newlevels <- colnames(newZAlist[[it]])
          if (identical(oldlevels,newlevels)) {
            oldnewC <- Diagonal(length(oldlevels)) ## replaces old identityMatrix
          } else {
            oldornew <- unique(c(oldlevels,newlevels))
            oldnewC <- diag(length(oldornew))
            colnames(oldnewC) <- rownames(oldnewC) <- oldornew
            oldnewC <- oldnewC[oldlevels,newlevels]
          }
          attr(oldnewC,"isEachNewInOld") <- newlevels %in% oldlevels  ## but this attr is unevaluated (-> NULL) for spatial models 
        }
        return(oldnewC)
      })
      #
      loclist <- list(X.pv=newX.pv,beta_w_cov=beta_w_cov,covMatrix=variances$cov,lambda=object$lambda)
      if (!is.null(uuCnewnew)) {
        uuCnewnewList <- lapply(seq_len(nrand),function(it) {
          if (it %in% spatialOne) { ## avoids running the next algo which is slow on large matrices
            return(uuCnewnew)
          } else return(NULL)
        }) ## should make sure it has the same structure as the other matrix lists
        loclist$CnewnewList <- uuCnewnewList 
      } ## else no loclist$CnewnewList (tested in calcPredVar) 
      #
      if ( ! is.null(newdata)) loclist$invColdoldList <- invColdoldList
      if (nrand==1L) {
        loclist$Coldnew <- oldnewClist[[1]]
        loclist$newZA <- newZAlist[[1]]
        ## but the code for nrand >1 should work for nrand==1:
      } else {
        loclist$Coldnew <- oldnewClist
        loclist$newZA <- newZAlist
      }
      if (variances$disp) loclist$logdispObject <- object$get_logdispObject()
      if (variances$cov) {
        respVar <- as.matrix(do.call(calcPredVar,loclist)) ## matrix, not Matrix (assumed below)
        rownames(respVar) <- colnames(respVar) <- rownames(locdata)
      } else {
        respVar <- do.call(calcPredVar,loclist) 
        names(respVar) <- rownames(locdata)
      }
    } else {
      if (variances$cov) {
        respVar <- matrix(0,nrow=nrow(locdata),ncol=nrow(locdata))
      } else respVar <- rep(0,nrow(locdata))
    }
  } else respVar <- rep(0,nrow(locdata))  
  if (! is.null(object$beta_cov)) {
    if ( variances$fixefVar || (nrand==0L && variances$linPred) ) {
      fixefcov <- newX.pv %*% object$beta_cov %*% t(newX.pv)
      if (variances$cov) {
        attr(resu,"fixefVar") <- fixefcov 
      } else attr(resu,"fixefVar") <- diag(fixefcov)
      if (nrand==0L) { ## otherwise there is already such a term in predVar
        respVar <- respVar + attr(resu,"fixefVar") 
      }
    }
  }
  attr(resu,"predVar") <- respVar ## vector or matrix
  if (variances$residVar) {
    if (object$family$family %in% c("poisson","binomial","COMPoisson")) {
      attr(resu,"residVar") <- object$family$variance(fv)
    } else attr(resu,"residVar") <- calcResidVar(object,newdata=locdata) 
    if (inherits(respVar,"matrix")) {
      diag(respVar) <- diag(respVar) + attr(resu,"residVar")
    } else respVar <- respVar + attr(resu,"residVar")
  }
  if (variances$respVar) attr(resu,"respVar") <- respVar
  if ( is.matrix(resu) && NCOL(resu)==1L) {
    class(resu) <- c("predictions",class(resu))
  } ## for print.predictions method which expects a 1-col matrix
  # intervals
  checkVar <- setdiff(intervals,names(attributes(resu)))
  if (length(checkVar)>0L) {
    warning(paste("Variances",paste(checkVar,collapse=", "),
                  "not available for interval computation.\n Check arguments."))
    intervals <- intersect(intervals,names(attributes(resu)))
  } 
  if(length(intervals)>0L) {
    intervalresu <- NULL
    for (st in intervals) {
      varcomp <- attr(resu,st)
      if (is.null(varcomp)) warning(paste("Prediction variance component",st,"requested but not available: check input."))
      if (is.matrix(varcomp)) varcomp <- diag(varcomp)
      eta <- object$family$linkfun(resu[,1L])
      pv <- 1-(1-level)/2
      sd <- qnorm(pv)*sqrt(varcomp)
      interval <- cbind(object$family$linkinv(eta-sd),object$family$linkinv(eta+sd))
      colnames(interval) <- paste(st,c(signif(1-pv,4),signif(pv,4)),sep="_")
      intervalresu <- cbind(intervalresu,interval)
    }
    attr(resu,"intervals") <- intervalresu
  }
  return(resu)
}

print.vcov.HLfit <-function(x, expanded=FALSE, ...) {
  a <- attributes(x)
  attr(x,"beta_v_cov") <- NULL  
  print.default(x)
  cat("with additional attribute(s):")
  std.attr <- c("names","dim","dimnames","class") ## attributes not to be shown
  nam <- names(a)
  if (expanded) { # shows structure of attributes as in utils:::str.default
    cat("\n")
    nest.lev <- 0
    indent.str <- paste(rep.int(" ", max(0, nest.lev + 1)), collapse = "..")
    strO <- getOption("str")
    strO <- modifyList(strOptions(), strO) ## seems to provide a format.fun function
    `%w/o%` <- function(x, y) x[is.na(match(x, y))]
    nfS <- names(fStr <- formals())
    ## this scans the substructure of each attribute
    strSub <- function(obj, ...) {
      nf <- nfS %w/o% c("object", "give.length", "comp.str", 
                        "no.list", names(match.call())[-(1:2)], "...")
      aList <- as.list(fStr)[nf]
      aList[] <- lapply(nf, function(n) eval(as.name(n)))
      strObj <- function(...) str(obj, ...)
      do.call(strObj, c(aList, list(...)), quote = TRUE)
    }
    for (i in seq_along(a)) if (all(nam[i] != std.attr)) {
      cat(indent.str, paste0("- attr(*, \"", nam[i], "\")="), 
          sep = "")
      strSub(a[[i]], give.length = TRUE, indent.str = paste(indent.str, 
                                                            ".."), nest.lev = nest.lev + 1)
    }
  } else {
    cat(" ")
    nam <- setdiff(nam,std.attr)
    cat(paste(nam,collapse=", "))  
    cat("\n")
  }
  invisible() ## do not return x since it has lost a useful attribute
}


`[.predictions` <- function (x, i, j, 
                             drop = TRUE ## by default, this function will return scalar/vector  
                             ) {
  class(x) <- "matrix" ## removes "predictions" => set back later
  #   if (is.data.frame(x)) {
  #     resu <- x[i,j]
  #   } else 
  resu <- x[i,j,drop=drop]
  if ( ! drop) {
    fixefVar <- attr(x, "fixefVar")
    if ( ! is.null(fixefVar)) {
      if (is.null(dim(fixefVar))) {
        fixefVar <- fixefVar[x]
      } else fixefVar <- fixefVar[x,x,drop=FALSE]
    }
    predVar <- attr(x, "predVar")
    if ( ! is.null(predVar)) {
      if (is.null(dim(predVar))) {
        predVar <- predVar[x]
      } else predVar <- predVar[x,x,drop=FALSE]
    }
    frame <- attr(x, "frame")
    if ( ! is.null(frame)) frame <- frame[x,] ## dataframe => nodrop
    residVar <- attr(x, "residVar")
    if ( ! is.null(frame)) residVar <- residVar[x,drop=FALSE]
    respVar <- attr(x, "respVar")
    if ( ! is.null(respVar)) {
      if (is.null(dim(respVar))) {
        respVar <- respVar[x]
      } else respVar <- respVar[x,x,drop=FALSE]
    }
    class(resu) <- c("predictions","matrix")
    structure(resu,fixefVar=fixefVar,predVar=predVar,residVar=residVar,frame=frame,fittedName=attr(x, "fittedName"))
  } else return(resu)
} # Use unlist() to remove attributes from the return value

print.predictions <- function (x, expanded=FALSE, ...) {
  asvec <- as.vector(x)
  rnames <- rownames(x)
  toolong <- nchar(rnames)>9
  rnames[toolong] <- paste(substr(rnames[toolong],0,8),".",sep="")
  names(asvec) <- rnames
  cat("Point predictions:\n")
  print(asvec)
  cat("*stored as* 1-col matrix with attributes:")
  std.attr <- c("names","dim","dimnames","class") ## attributes not to be shown
  a <- attributes(x)
  nam <- names(a)
  if (expanded) { # shows structure of attributes as in utils:::str.default
    cat("\n")
    nest.lev <- 0
    indent.str <- paste(rep.int(" ", max(0, nest.lev + 1)), collapse = "..")
    strO <- getOption("str")
    strO <- modifyList(strOptions(), strO) ## seems to provide a format.fun function
    `%w/o%` <- function(x, y) x[is.na(match(x, y))]
    nfS <- names(fStr <- formals())
    ## this scans the substructure of each attribute
    strSub <- function(obj, ...) {
      nf <- nfS %w/o% c("object", "give.length", "comp.str", 
                        "no.list", names(match.call())[-(1:2)], "...")
      aList <- as.list(fStr)[nf]
      aList[] <- lapply(nf, function(n) eval(as.name(n)))
      strObj <- function(...) str(obj, ...)
      do.call(strObj, c(aList, list(...)), quote = TRUE)
    }
    for (i in seq_along(a)) if (all(nam[i] != std.attr)) {
      cat(indent.str, paste0("- attr(*, \"", nam[i], "\")="), 
          sep = "")
      strSub(a[[i]], give.length = TRUE, indent.str = paste(indent.str, 
                                                            ".."), nest.lev = nest.lev + 1)
    }
  } else {
    cat(" ")
    nam <- setdiff(nam,std.attr)
    cat(paste(nam,collapse=", "))  
    cat("\n")
  }
  invisible()
}