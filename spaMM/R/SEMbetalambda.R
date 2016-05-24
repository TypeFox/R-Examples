## y* = X beta + Z b + e
## b :=autocorrelated ranefs, =Lv for iid v, = Uw for i but not id w in CAR
SEMbetalambda <- function(beta_eta,lambda,
                          corr_est=NULL,
                          lambda.Fix=NULL, 
                          nSEMiter=200, # must be >9 cf check of non-default value in preprocess
                          ngibbs=20,
                          nMCint=10000, ## 10000 see notes from 11/05/2015
                          SEMseed=NULL,SEMsample=NULL,SEMlogL="pmvnorm",
                          control_pmvnorm=list(),  ## 'for development purposes, not documented...'
                          control_pMVN=list(),  ## 'for development purposes, not documented...'
                          ZA,X.pv,qr.XtX,
                          symSVD, ## 
                          ZAL,whichy1,off,
                          stop.on.error,
                          verbose=FALSE,
                          X_lamres,
                          control.glm,
                          mc){ 
  if (FALSE) {
    adaptive <- TRUE
    nSEMiter <- 10000L
    pilotLength <- 50L
  } else {
    adaptive <- FALSE
    pilotLength <- nSEMiter
  }
  if (is.null(SEMsample)) SEMsample <- ceiling(nSEMiter/2):nSEMiter
  if(!is.null(SEMseed)) { ## direct HLCor call
    set.seed(SEMseed) ## so that estimates of beta,lambda are repeatable ## comment ne pas avoir a retirer les memes nombres XXXX fois ?
    #      cat(paste("SEMseed=",SEMseed))
  } # else HLCorobj call
  nobs <- nrow(X.pv)
  betaMat <- matrix(0,nrow=nSEMiter,ncol=ncol(X.pv))
  colnames(betaMat) <- colnames(X.pv) 
  betaMat[1,] <- beta_eta
  lambdaVec <- numeric(nSEMiter)
  condVar <- rep(0,nSEMiter)
  condVar[1] <- lambdaVec[1] <- lambda
  corr.model <- symSVD$corr.model  
  estimCARrho <- ( corr.model=="adjacency" &&  ! is.null(corr_est$rho)) 
  if (estimCARrho) {
    rhoVec <- numeric(nSEMiter)
    rhoVec[1] <- corr_est$rho
  }
  n_u_h <- symSVD$dim[2L] ## locally true as implied by code below.
  EbGivenY <- matrix(0,nrow=nSEMiter,ncol=n_u_h)
  decomp <- symSVD$symsvd
  if ( ! is.identity(ZA)) {
    ZAisI <- FALSE
    ZA <- as.matrix(ZA) ## useless to stay Matrix
    ZAE <- ZA %id*id% decomp$u ## dense in spatial model  
  } else ZAisI <- TRUE
  whichy0 <- (! whichy1) ##FR->FR in preprocess ?
  ny1 <- sum(whichy1) ##FR->FR in preprocess ?
  ny0 <- nobs - ny1
  if (SEMlogL == "pmvnorm") {
    pmvnorm.lower <- rep(0,nobs)
    pmvnorm.lower[whichy0] <- -Inf
    pmvnorm.upper <- rep(0,nobs)
    pmvnorm.upper[whichy1] <- Inf
  }
  gibbsSample <- ceiling(ngibbs/2):ngibbs
  # prevmsglength <- 0L
  i <- 2
  while (TRUE) {
    ## whatever depends on rhoVec[i-1] (fixed in the gibbs block)
    if ( estimCARrho ) {
      decomp$d <- 1/(1-rhoVec[i-1]*decomp$adjd) ## le decomp$d originel est Corr, pas pour L
    }
    if (i==2L || estimCARrho) { 
      if (corr.model=="identity") {
        LMatrix <- invLMatrix <- decomp$u
      } else {
        LMatrix <- ZWZt(decomp$u,sqrt(decomp$d))
        invLMatrix <- ZWZt(decomp$u,1/sqrt(decomp$d)) ## no wrapper (testing identity is a loss of time)
      }
      if( ! ZAisI) {
        ZAEdEAZ <- ZWZt(ZAE,decomp$d) ## FR->FR large (nresp * nresp)
        forV <- selfAdjointSolverCpp(ZAEdEAZ) ## so that inv(V) = ZWZt(forV$u,1/(1+lambda_est * forV$d))
        ##             [ D = lambda Corr]  . Z'      . forV$u but without the lambda factor
        if (corr.model=="identity") {
          ranefCorr <- decomp$u
        } else {
          ranefCorr <- tcrossprodCpp(LMatrix) ## reconstruction compliqu?e de la corr matrix ?
        }
        LHSCorrblob <- as.matrix(ranefCorr %id*% t(ZA) %*% forV$u)
      }
    }
    ## whatever (also) depends on lambdaVec[i-1] (fixed in the gibbs block)
    if(ZAisI) {
      CondNorm <- CondNormfn(decomp,lambdaVec[i-1])
    } else { ## D - D Z' inv(V) Z D
      ##          [D = lambda *Corr] - [LHSblob= lambda Corr Z' forV$u]. 1/(1+lambda_est * forV$d) . t(LHSblob) 
      ## with lambda^2 /(1+lambda d) = lambda/(1/lambda + d)
      condCov <- lambdaVec[i-1] * ranefCorr - ZWZt(LHSCorrblob,lambdaVec[i-1]/(1/lambdaVec[i-1] + forV$d)) ## ZWZt must be a slow step
      condL <- RcppChol(as.matrix(condCov))$L ## such that only tcrossprod(condL) = tcrossprod(tcrossprod(condL)) when ZAisI
      ## not more code because I will try to perform only matrix * vector operations
    }
    # S part of the SEM algorithm
    # we use a Gibbs sampling algorithm
    randbGivenObs <- sqrt(lambdaVec[i-1]) * (LMatrix %*% rnorm(n_u_h,0))
    augY <- rep(0,nrow(ZAL))
    fix <- X.pv %*% betaMat[i-1,] + off
    if (estimCARrho) {
      condw2 <- rep(0,n_u_h)      
    } else condv2 <- rep(0,n_u_h)
    Estcond_bMeans <- rep(0,n_u_h) ## denoting the fact that the conditioning event is itself random in the gibbs
    for (k in 1:ngibbs) {
      # random generation of augY given obs: y and v (fixed beta, fixed lambda)
      moy.augY <- fix + ZA %*% randbGivenObs  
      augY[whichy1] <- rntpos(ny1,moy.augY[whichy1],1)
      augY[whichy0] <- rntneg(ny0,moy.augY[whichy0],1)
      ## whatever depends on augY
      if(ZAisI) {
        randcond_bMean <- CondNorm$condLvReg %*% (augY-fix)
        randcond_brand <- CondNorm$sqrtCondCovLv %*% rnorm(n_u_h,0)
      } else { ##randcond_bMean <- lambda_est * (ranefCorr %*% (t(ZA) %*% solve(augYCov,augY-fix))) ## DZ'inv(V)(y-X beta) in Searle p. 275
        ##             [D = lambda *Corr]   .     Z'    .   inv(V).(augY-fix)   with initial lambda brought inside
        ##randcond_bMean <- ranefCorr %*% t(t(forV$u %*% t((t(augY-fix) %*% forV$u)/(1/lambdaVec[i-1] + forV$d))) %*% ZA) ## only t(vector) ## SLOW
        ##randcond_bMean <- LHSCorrblob %*% t((t(augY-fix) %*% forV$u)/(1/lambdaVec[i-1] + forV$d)) ## only t(vector)
        # faste code for the same:
        locv <- augY-fix
        dim(locv) <- c(1,nrow(locv)) ## fast transposition
        locv <- (locv %*% forV$u)/(1/lambdaVec[i-1] + forV$d)
        dim(locv) <- c(ncol(locv),1) ## fast transposition
        randcond_bMean <- LHSCorrblob %*% locv
        randcond_brand <- condL %*% rnorm(n_u_h,0) ## has zero expectation
      }
      ## augY should be fix + randbGivenObs + one-epsilon-per-individual 
      # random generation of v given (y and) augmented Y 
      randbGivenObs <- randcond_bMean + randcond_brand ## b
      if (k %in% gibbsSample) {
        Estcond_bMeans <- Estcond_bMeans + randcond_bMean ## sufficient for E[b|obs] since randcond_brand has zero expectation
        if (estimCARrho) {
          condw2 <- condw2 + (t(decomp$u) %*% randbGivenObs)^2 ## w indep mais pas i.d., pour estim lambda et rho
        } else condv2 <-  condv2 + (invLMatrix %*% randbGivenObs)^2 ## v iid
      }
    } ## end ngibbs loop
    ### end of E step
    ### M step of the SEM algorithm
    # determination of beta by standard least square #betaMat[i,] <- lm((z-Lvs)~X.pv-1)$coeff
    EbGivenY[i,] <- Estcond_bMeans/length(gibbsSample) ##  E[b|...]
    if (ZAisI) {
      betaMat[i,] <- solveWrap.vector( qr.XtX , t(X.pv) %*% (augY- EbGivenY[i,] -off) ,stop.on.error=stop.on.error)
    } else betaMat[i,] <- solveWrap.vector( qr.XtX , t(X.pv) %*% (augY-(ZA %*% EbGivenY[i,]) -off) ,stop.on.error=stop.on.error)
    # Estim of (co)variance parameters
    if (is.null(lambda.Fix)) {
      if ( estimCARrho ) {
        condw2 <- condw2/length(gibbsSample)
        locdf <- data.frame(v2=condw2,adjd=decomp$adjd) ## FR->FR hmf. adjd as part of an X_lamres computed by preprocess ?
        ## these two functions are (equivalently) maximized by the GLM:
        #         bordel <- function(lr) {
        #           covmat <- diag(lr[1]/(1-lr[2]*decomp$adjd))
        #           dmvnorm(as.numeric(vs), mean = rep(0, nrow(covmat)), sigma = covmat, log = TRUE)
        #         }
        #         bordel2 <- function(lr) {
        #           covmat <- lr[1] * ZWZt(decomp$u,1/(1-lr[2]*decomp$adjd)) ## pour randbGivenObs
        #           dmvnorm(as.numeric(randbGivenObs), mean = rep(0, nrow(covmat)), sigma = covmat, log = TRUE)
        #         }        
        #print(fixef(HLfit(v2 ~ 1+ adjd,family=spaMM_Gamma("inverse"),prior.weights=rep(1/2,n_u_h),data=locdf)))
        resglm <- suppressWarnings(glm(v2 ~ 1+ adjd,family=spaMM_Gamma("inverse"),
                                       weights=rep(1/2,n_u_h),data=locdf,start=c(1/mean(locdf$v2),0))) ## is ML (no EQL issue)
        coeffs <- coefficients(resglm) 
        lambdaVec[i] <- 1/coeffs[1]
        rhoVec[i] <- - coeffs[2]/coeffs[1] ## note that lambdas[k]= mean((1-rhos[k]*decomp$adjd)*vs^2)
        if (verbose) prevmsglength <- overcat(paste("iter", i,": lambda=",signif(lambdaVec[i],4),"; rho=",signif(rhoVec[i],4)), prevmsglength)
      } else {
        condv2 <- condv2/length(gibbsSample)
        lambdaVec[i] <- mean(condv2)
        if (verbose) prevmsglength <- overcat(paste("iter",i,": lambda=",signif(lambdaVec[i],4)), prevmsglength)
      }
    } else {
      if ( estimCARrho ) {
        condw2 <- condw2/length(gibbsSample)
        locdf <- data.frame(v2=condw2,adjd=decomp$adjd) ## FR->FR hmf. adjd as part of an X_lamres computed by preprocess ?
        resglm <- suppressWarnings(glm(v2 ~ adjd -1,offset=rep(1/lambda.Fix,nrow(locdf)),family=spaMM_Gamma("inverse"),
                                       weights=rep(1/2,n_u_h),data=locdf,start=c(0))) ## is ML (no EQL issue)
        coeffs <- coefficients(resglm) 
        lambdaVec[i] <- lambda.Fix
        rhoVec[i] <- - coeffs[1]*lambda.Fix 
        if (verbose) prevmsglength <- overcat(paste("iter", i,": rho=",signif(rhoVec[i],4)), prevmsglength)
      } else lambdaVec[i] <- lambda.Fix
    }
    if (adaptive && i==pilotLength) {
      ## see memo in adaptiveDiagnostics.R.txt
    } 
    if (nSEMiter==i) { ## either !adaptive or (adaptive && nSEMiter=pilotLength=i)
      break;
    } #else { 
      #if (nSEMiter < i) {
      #  lambdaVec <- lambdaVec[seq_len(i)] ## reduction
      #} else {
      #  oldvec <- lambdaVec
      #  lambdaVec <- numeric(nSEMiter)
      #  lambdaVec[seq_len(i)] <- oldvec
      #}
    #} 
    i <- i+1
  }  ## end nSEMiter loop
  if (verbose) {
    cat("\n") 
    plot(lambdaVec[1:nSEMiter])
    if (estimCARrho) plot(rhoVec)
  }
  
  logL_from_SEMsample <- function(SEMsample) {
    beta_eta <- colMeans(betaMat[SEMsample,,drop=FALSE]) 
    lambda <- mean(lambdaVec[SEMsample]) ## lambda_est will be given a different length
    if (estimCARrho) {
      corr_est$rho <- mean(rhoVec[SEMsample])
      decomp$d <- 1/(1-corr_est$rho*decomp$adjd) 
      LMatrix <- ZWZt(decomp$u,sqrt(decomp$d))
      invLMatrix <- ZWZt(decomp$u,1/sqrt(decomp$d))
      if (ZAisI) {ZAL <- LMatrix} else { ZAL <- ZA %id*id% LMatrix }
    }
    ## estim lik
    fix <- X.pv %*% beta_eta + off
    if (verbose) cat("Estimating the likelihood...")
    method <- SEMlogL
    if (SEMlogL %in% c("halton","pMVN","pmvnorm")) { 
      ZALtZAL <- tcrossprodCpp(ZAL)
      pmvnorm.Sig <- lambda*ZALtZAL+diag(nobs)
    }
    if (SEMlogL=="pmvnorm") { 
      pmvnorm_arglist <- c(control_pmvnorm,
                           list(lower=pmvnorm.lower,upper=pmvnorm.upper,mean=as.vector(fix),sigma=pmvnorm.Sig,
                                # mvt.f shows that the algo continues (within some limits) 
                                #  if the error is higher than the MAX of the abseps and releps-induced constraints. 
                                # For large negative logL, we can only control releps, and we must set set abseps=0 to control releps.
                                abseps=0,releps=0.5))
      if (is.null(pmvnorm_arglist$maxpts)) ##  postprocessing in corrHLfit provides a maxpts
        pmvnorm_arglist$maxpts <- quote(1000L*nrow(pmvnorm.Sig)) ## "sensible ... to start with MAXPTS = 1000*N" (mvt.f)
      ## but the default in optim through smooth is 250L*nrow(pmvnorm.Sig)
      # CALL
      Lapp <- do.call("pmvnorm",pmvnorm_arglist) 
      # 
      # pmvnorm may fail in various ways. Worst,
      ## pmvnorm bug, difficult to replicate without control of seed:
      ## Lapp my be NaN with "normal completion" message.
      pmvnorm_fatal <- (is.na(Lapp) || Lapp[1]==0 || attr(Lapp,"error")==0) 
      if ( (! pmvnorm_fatal) && attr(Lapp,"msg")=="Completion with error > abseps" ) {
        ## we could increase $maxpts, but this does not appear useful (and must handle the quote object)
        # we could now use $abseps > 0 but given the $releps=0.5 failed, it's not clear what this would improve
        # => do nothing, keep the large Lapp error
      }
    } else pmvnorm_fatal <- FALSE
    if (pmvnorm_fatal || SEMlogL %in% c("halton","pMVN")) { ## input rand_seq or default GHK method
      Lmat <- RcppChol(pmvnorm.Sig)$L ## t(R::chol)
      pMVN_arglist <- c(control_pMVN, 
                    list(L=Lmat,limits = as.vector(fix),ismax=whichy1))
      #if (SEMlogL=="halton") { ## but avoid the fatal qrng::ghalton() for dim > 360
      #  pMVN_arglist$rand_seq <- randtoolbox::halton(4L*nrow(pmvnorm.Sig),d=nrow(pmvnorm.Sig),init=FALSE,usetime=FALSE)
      #  pMVN_arglist$nrep <- NULL
      #} else 
      if (is.null(pMVN_arglist$nrep)) 
        pMVN_arglist$nrep <- quote(200L*nrow(pmvnorm.Sig)) ## cf notes 22/03/2016
      # CALL
      blob <- do.call("pMVN",pMVN_arglist) ## GHK by default
      # 
      logLapp <- blob$logInt
      seInt <- blob$seInt
      attr(logLapp,"method") <- "      logL (GHK)" ## directly usable for screen output
    } else if (SEMlogL=="p_v") { ## use standard Laplace approx for estimating the likelihood
      ## this appears to perform very poorly probably due to poor estim of lambda by Laplace estimation.
      ## HACK for performing the computation: extracts any 'processed' argument from the call, substitutes PQL/L to SEM in it, etc. 
      arglist <- as.list(mc) ## of which [[1]] is "HLfit"
      proc <- arglist$processed 
      if (! is.null(proc)) {
        arglist$processed <- eval.update.call(proc$callargs,HLmethod="PQL/L") ## update 'preprocess' call
        arglist$HLmethod <- NULL ## for clarity; should be ignored anyway as the processed$HL info should be used
      } else arglist$HLmethod <- "PQL/L"
      arglist$etaFix$beta <- beta_eta
      #      arglist$etaFix$v_h <- v_h ## interestingly, disastrous. Probably also related to poor estim of lambda by Laplace approx
      arglist$ranFix$lambda <- lambda
      if (estimCARrho) {
        stop("code needed to (effectively) update ZAL in arglist ...")
      }
      logLapp <- eval(as.call(arglist))$APHLs$p_v
      seInt <- 0
      attr(logLapp,"method") <- "p_v(h) (marginal L):"
    } else if (SEMlogL=="pmvnorm") {
      logLapp <- log(Lapp[1L])
      seInt <- attr(Lapp,"error")/(2.575829*Lapp) ## = qnorm(0.995) because pmvnorm's error is def'ine'd as bound of 99% interval
      attr(logLapp,"method") <- "  logL (pmvnorm)" ## directly usable for screen output
    } else stop(paste("Unknown integration method",SEMlogL,"requested"))
    ## naive MC removed 04/2016, GHK is now reliable
    return(list(beta_eta=beta_eta,lambda=lambda,corr_est=corr_est,logLapp=logLapp,seInt=seInt))
  }
  
  resu <- logL_from_SEMsample(SEMsample)
  lambda <- resu$lambda ## snif... 2015/09/16
  
  #   #jackknife (inefficient procedure recomputing matrix decomp...)
  #   jackmean <- sapply(seq_len(length(SEMsample)),function(it) {logL_from_SEMsample(SEMsample[-it])$logLapp})
  #   meanjackmean <- mean(jackmean)
  #   nextlogLapp <- resu$logLapp + (length(SEMsample)-1L)*(resu$logLapp-meanjackmean) 
  #   resu$logLapp <- nextlogLapp
  #   # totally inefficient du to variance in estimation...
  
  if (verbose) cat(" SEM done.\n")
  if (estimCARrho) {
    ## to provide a conforming glm to HLfit but also coefficients to be substituted to the coefficents of this glm
    locdf$resp <- condw2
    resu$glm_lambda <- calc_CARdispGammaGLM(data=locdf, lambda.Fix=lambda.Fix,control=control.glm)
    resu$glm_lambda$coeffs_substitute <- c("(Intercept)"= 1/lambda,"adjd"= - resu$corr_est[["rho"]]/lambda)
  } else {
    ## calc_dispGammaGLM needs 'data' even if 'data' has zero cols  
    Intcol <- which(colnames(X_lamres)=="(Intercept)")
    loclamdata <- data.frame(X_lamres[, - Intcol,drop=FALSE])
    ## FR->FR in a "later" version the formula and 'data' should be provided by processed$ as X_lamres come from there
    ## FR->FR could be tested now ?
    if (ncol(loclamdata)>0L) {
      loclamformula <- as.formula(paste("~",paste(colnames(loclamdata),collapse="+")))
    } else loclamformula <- ~1
    resu$glm_lambda <- calc_dispGammaGLM(
      formula=loclamformula,
      dev.res=as.matrix(condv2),
      lev=rep(0,n_u_h),
      data=loclamdata,
      etastart=rep(log(lambda),n_u_h), 
      control=control.glm)
    resu$glm_lambda$coeffs_substitute <- c("(Intercept)"= log(lambda))
  }
  ## output for prediction only: 
  resu$v_h <- invLMatrix %*% colMeans(EbGivenY[SEMsample,,drop=FALSE]) ## v iid from b not iid
  return(resu)
} ## end def of SEMbetalambda

