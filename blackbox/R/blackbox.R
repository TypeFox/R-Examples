bboptim <-function(data,ParameterNames=NULL,respName=NULL,control=list(),force=FALSE,optimizers=blackbox.getOption("optimizers")) {
  np <- ncol(data)-1L
  if (is.null(ParameterNames)) ParameterNames <- colnames(data)[seq(np)]
  if (is.null(respName)) respName <- colnames(data)[np+1L]
  sorted_etc <- prepareData(data=data,ParameterNames=ParameterNames,
                      respName=respName,verbose=FALSE) ## sets min- and maxSmoothness
  #
  fnscale <- control$fnscale
  maximizeBool <- control$maximize
  if (is.null(maximizeBool)) {
    maximizeBool <-  ! ( is.null(fnscale) || fnscale >= 0 ) ## default is FALSE
  } else if (is.null(fnscale)) {
    if (maximizeBool) fnscale <- -1
  } else stop("Mixing controls maximize and fnscale is dangerous. Please correct.")
  gcvres <- calcGCV(sorted_etc, force=force, decreasing=maximizeBool, optimizers=optimizers ) ## decreasing important info for calcGCV -> selectFn
  ranfix <- list(rho=1/gcvres$CovFnParam[ParameterNames],
                 nu=gcvres$CovFnParam[["smoothness"]])
  form <- as.formula(paste(respName,"~ 1 + Matern(1|",paste(ParameterNames,collapse=" + "),")"))
  optimizers[optimizers=="optim"] <- "L-BFGS-B"
  # Re-estimation of lambda and phi, cf comments in Infusion
  # Krig_coef is not called, hence no new CSmooth object is allocated, hence there's nothing to delete later
  phi <- max(gcvres$pureRMSE^2,1e-8) ## FR->FR spaMM accepts low phi values, corrects them only in cal.p_v
                                      # this is confusing as clik is then not stable to translation...
  ranfix <- c(ranfix,list(phi=phi))
  if (TRUE) {
    thisfit <- HLCor(form,data=data,ranPars=ranfix,REMLformula=form)
    ranfix <- c(ranfix,list(lambda=thisfit$lambda))
    etafix <- list(beta=fixef(thisfit))
  } else { ## FR->FR equivalent alternatives? except that one provides etaFix the other not
    ranfix <- c(ranfix,list(lambda=phi/gcvres$lambdaEst)) ## lambdaGCV := phiHGLM/lambdaHGLM !
    etafix <- list()
  }
  thisfit <- corrHLfit(form,data=sorted_etc,ranFix=ranfix,etaFix=etafix) ## full data smoothed with lambda and phi estimated from cleaned data
  #
  lower <- apply(sorted_etc[,ParameterNames,drop=FALSE],2,min)
  upper <- apply(sorted_etc[,ParameterNames,drop=FALSE],2,max)
  if ( ! maximizeBool ) {
    init <- thisfit$data[which.min(predict(thisfit)),ParameterNames]
    optr_fitted <- list(par=init,value=min(predict(thisfit)))
  } else {
    init <- thisfit$data[which.max(predict(thisfit)),ParameterNames]
    optr_fitted <- list(par=init,value=max(predict(thisfit)))
  }
  mse <- attr(predict(thisfit,init,variances=list(linPred=TRUE,dispVar=TRUE)),"predVar")[1]
  mse <- max(0,mse)
  optr_fitted$RMSE <- sqrt(mse)
  ## in all cases, try to keep the structure of the return object for each method, but add easily accessible members
  if ("L-BFGS-B" %in% optimizers) {
    if (is.null(control$parscale)) control$parscale <- (upper-lower)
    optr <- optim(init, function(v) {as.numeric(predict(thisfit,newdata=v))},
                  ## as numeric because otherwise in 1D, optim -> minimize -> returns a max
                  ##   of same type as predict(thisfit$fit,newdata=v), i.e. matrix.
                  lower=lower,upper=upper,control=control,method="L-BFGS-B")
  } else if ("lbfgsb3" %in% optimizers){
    if ( maximizeBool ) {
      ufn <- function (v) { - as.numeric(predict(thisfit,newdata=v)) }
    } else {
      ufn <- function (v) { as.numeric(predict(thisfit,newdata=v)) }
    }
    if ( ! requireNamespace("lbfgsb3",quietly=TRUE) ) {
      stop("Package lbfgsb3 not installed.")
    }
    optr <- lbfgsb3::lbfgsb3(prm=unlist(init),fn=ufn,lower=lower,upper=upper)
    if ( maximizeBool ) {optr$value <- - optr$f} else {optr$value <- optr$f}
    optr$f <- NULL
    optr$par <- optr$prm; optr$prm <- NULL
    #optr$counts <- optr$info$isave[34]
  } else if ("bobyqa"  %in% optimizers) {
    if ( maximizeBool ) {
      ufn <- function (v) { - as.numeric(predict(thisfit,newdata=v)) }
    } else {
      ufn <- function (v) { as.numeric(predict(thisfit,newdata=v)) }
    }
    control <- list(rhobeg=min(abs(upper-lower))/20)
    control$rhoend <- max(1,control$rhobeg)/1e6
    if ( ! requireNamespace("minqa",quietly=TRUE) ) {
      stop("Package minqa not installed.")
    }
    optr <- minqa::bobyqa(par=unlist(init),fn=ufn,lower=lower,upper=upper,control=control)
    if ( maximizeBool ) {optr$value <- - optr$fval} else {optr$value <- optr$fval}
  } else { ## default
    if ( maximizeBool ) {
      ufn <- function (v) { - as.numeric(predict(thisfit,newdata=v)) }
    } else {
      ufn <- function (v) { as.numeric(predict(thisfit,newdata=v)) }
    }
    optr <- nloptr(x0=unlist(init),eval_f=ufn,lb=lower,ub=upper,
                   opts=list(algorithm="NLOPT_LN_BOBYQA",maxeval=-1))
    if ( maximizeBool ) {optr$value <- - optr$objective} else {optr$value <- optr$objective}
    optr$par <- optr$solution
  }
  colTypes <- list(ParameterNames=ParameterNames,respName=respName)
  eta <- predict(thisfit,newdata=optr$par,variances=list(linPred=TRUE,dispVar=TRUE))
  RMSE <- sqrt(max(0,attr(eta,"predVar")))
  # assessment of convergence
  reltol <- control$reltol
  if (is.null(reltol)) reltol <- sqrt(.Machine$double.eps)
  if (maximizeBool) {
    convergence <- (optr_fitted$value > optr$value-reltol)
  } else convergence <- (optr_fitted$value < optr$value+reltol)
  resu <- list(fit=thisfit,optr_fitted=optr_fitted,optr=optr,callArgs = as.list(match.call())[-1],colTypes=colTypes,
               GCVmethod=gcvres$GCVmethod,RMSE=RMSE,convergence=convergence)
  class(resu) <- c("list","bboptim")
  return(resu)
}

print.bboptim <- function(x,...) {
  if (x$convergence) {
    cat("\nApparent convergence at optimum:\n")
    print(x["optr"],...)
    cat("Other elements of bboptim object not fully shown:",
        paste(setdiff(names(x),c("optr","convergence")),collapse=", "))
  } else {
    cat("\nApparent optimum not ascertained:\n")
    cat("Inferred optimum among fitted points:\n")
    print(x[["optr_fitted"]],...)
    cat("Inferred optimum over all space:\n")
    print(c(x[["optr"]][c("par","value")],RMSE=signif(x$RMSE,4)),...)
    cat("Other elements of bboptim object not fully shown:",
        paste(setdiff(names(x),c("optr_fitted","RMSE","convergence")),collapse=", "))
  }
}

summary.bboptim <- function(object,...) {
  print(x=object,...)
}

bbrhull <- function(n, from=10*n, vT,
                    object, ## used for predictor and for Qmax/Qmin
                    fixed=NULL, outputVars=object$colTypes$ParameterNames) {
  trypoints <- data.frame(rvolTriangulation(from, vT))
  colnames(trypoints) <- colnames(vT$vertices) ## supposeque non null...
  if (! is.null(fixed)) {
    trypoints <- cbind(trypoints, fixed)
  }
  if ("Migraine" %in% blackbox.getOption("usedBy")) {
    trypoints <- apply(trypoints, 1, tofullKrigingspace)
    if (length(outputVars)>1L) {
      trypoints <- t(trypoints)
    } else trypoints <- matrix(trypoints,ncol=1L)
  }
  colnames(trypoints) <- outputVars ## 'apply' feature
  trypred <- predict(object, newdata=as.data.frame(trypoints), predVar=TRUE)
  trySE <- attr(trypred, "predVar")
  trySE[trySE<0] <- 0
  trySE <- sqrt(trySE)
  Qmax <- object$Qmax
  if (is.null(Qmax)) {## minimization
    Qmin <- object$Qmin
    tryQ <- trypred - 1.96*trySE ## for candidate points
    expectedImprovement <- trySE*dnorm((tryQ-Qmin)/trySE)+(Qmin-tryQ)*pnorm((Qmin-tryQ)/trySE) ## 7.5 p. 121
  } else { ## maximization
    ## note difference with migraine application which switches from - logL to logL here
    tryQ <- trypred + 1.96*trySE ## for candidate points
    expectedImprovement <- trySE*dnorm((Qmax-tryQ)/trySE)+(tryQ-Qmax)*pnorm((tryQ-Qmax)/trySE) ## 7.5 p. 121
  }
  trypoints <- trypoints[order(expectedImprovement, decreasing=TRUE)[seq_len(n)], outputVars, drop=FALSE]
  return(trypoints) ## with names of the restricted space represented by vT
}

rbb <- function(object,n=NULL,from=NULL,focus=0.75) {
  if (focus<0 || focus>1) {
    focus <- max(0,min(1,focus))
    message.redef(paste("'focus' outside (0,1), set to",focus))
  }
  ParameterNames <- object$colTypes$ParameterNames
  np <- length(ParameterNames)
  if (is.null(n)) {
    n <- floor(10*(1+3*log(np))) ## 10 30 42 51 58 63 68 72 75 79
    n <- min(n,2^(np+1))         ##  4  8 16 32 58 63 ... ## rbb called only for resamples, not for initial grid
  }
  if (is.null(from)) {
    from <- n* floor(10*(1+3*log(np))) ## n* 10 30 42 51 58 63 68 72 75 79
    ## ## full default was 40 240 672 ... 6241 for np=10
  }
  if (from<=2*n) stop("from < 2*n: please increase 'from' relative to 'n'")
  ycolname <- object$colTypes$respName
  obspred <- predict(object$fit,variances=list(linPred=TRUE,dispVar=TRUE),binding=ycolname)
  obsSE <- attr(obspred,"predVar")
  obsSE[obsSE<0] <- 0 ## anticipating numerical problems
  fnscale <- object$callArgs$control$fnscale
  if ( ! is.null(fnscale) && fnscale < 0) { ##  maximization
    object$fit$Qmax <- max(obspred[,ycolname]+1.96 * sqrt(obsSE)) ## best improvement for already computed points
    ordr <- order(obspred[,ycolname],decreasing=TRUE)
  } else { ## minimization
    object$fit$Qmin <- min(obspred[,ycolname]-1.96 * sqrt(obsSE)) ## best improvement for already computed points
    ordr <- order(obspred[,ycolname],decreasing=FALSE)
  }
  # sampling from the whole parameter space (using object$fit$Qmin/max)
  vT <- volTriangulationWrapper(object$fit$data[,ParameterNames,drop=FALSE])
  ## used from=from/2 up to version 0.4.31:
  candidates <- bbrhull(round(n*(1-focus)), from=from,vT, object$fit,outputVars=ParameterNames) ## ordered points
  #
  if (FALSE) {## used up to version 0.4.31
    ntop <- min(nrow(obspred),4*np)
    vT <- volTriangulationWrapper(object$fit$data[ordr[seq(ntop)],ParameterNames,drop=FALSE])
    candidates2 <- bbrhull(n-nrow(candidates)-1L, from=from/2,vT, object$fit,outputVars=ParameterNames) ## ordered points
  } else { ## from version 0.4.32, 2015/11/17
    whichSimplex <- locatePointinvT(object$optr$par,vT)
    if (length(whichSimplex)==1L) {
      ## more transparent, but the general 'else' code applies...
      simplex <- vT$vertices[vT$simplicesTable[whichSimplex,],,drop=FALSE]
      candidates2 <- replicate(n-nrow(candidates)-1L,
                                 rsimplex(simplex=simplex))
      if (np == 1L) {
        candidates2 <- as.matrix(candidates2)
      } else candidates2 <- t(candidates2)
      rownames(candidates2) <- NULL
      colnames(candidates2) <- ParameterNames ## not NULL anyway
    } else {
      subvT <- subsimplices.volTriangulation(vT,whichSimplex)
      candidates2 <- rvolTriangulation(n=n-nrow(candidates)-1L,subvT)
    }
  }
  candidates <- rbind(candidates,candidates2,object$optr$par)
  # avoiding more than two replicates
  candidatesWPred <- predict(object$fit,newdata=candidates,binding=ycolname)
  poolWresp <- rbind(object$fit$data,candidatesWPred) ## mixes observed responses and predicted responses
  uli <- ULI(poolWresp[,1:np,drop=FALSE])
  ## test uli without optr$par
  if (any(table(uli[-length(uli)])>2L)) { ## on previous + candidates: should never be reached
    if (blackbox.getOption("dump_frames")) { ## FR->FR private option for debugging, removable later ?
    warning("Some candidate(s) in >2 copies.")
    dump.frames(dumpto = "dump_in_rbb",to.file=TRUE)
    } else stop("Some candidate(s) in >2 copies.")
  }
  ## tests with optr$par
  table_uli <- table(uli)
  curr_nrep_optr <- table_uli[uli[length(uli)]]
  if (curr_nrep_optr==1L) { ## add one more optr$par copy
    candidates <- rbind(candidates,object$optr$par)
    ndoublons <- 0L
  } else if (curr_nrep_optr==3L) { ## remove the latest optr$par copy and replicates two other point
    candidates <- candidates[-nrow(candidates),,drop=FALSE]
    ndoublons <- 2L
  } else if (curr_nrep_optr>3L) { ## this should never be reached
    stop("Candidate optimum already in >2 copies in candidate points.")
  } else ndoublons <- 1L ## curr_nrep_optr==2
  if (ndoublons>0L) {
    ## add a copy of the best predicted or observed non-replicated point
    singletonsWresp <- poolWresp[which(uli %in% which(table_uli==1L)),,drop=FALSE]
    if (nrow(singletonsWresp)>0L) {
      whichPts <-  order(singletonsWresp[,ycolname],
                         decreasing=( ! is.null(fnscale) && fnscale < 0))[seq_len(ndoublons)]
      candidates <- rbind(candidates,singletonsWresp[whichPts,1:np])
    } ## else odd but possible case of no singletons, ignore...
  }
  if ( ! inherits(candidates,"data.frame")) {
    if (blackbox.getOption("dump_frames")) { ## FR->FR private option
    warning("! inherits(candidates,\"data.frame\").")
    dump.frames(dumpto = "another_dump_in_rbb",to.file=TRUE) ## FR->FR private option for debugging, removable later ?
    } else stop("! inherits(candidates,\"data.frame\").")
  }
  candidates
}

