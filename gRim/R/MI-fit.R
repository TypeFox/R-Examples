###
### fit method for mModel objects
###
fit.mModel <- function(object, method="general", details=0, eps.parm=1e-6, maxit=100,...){
  method <- match.arg(method, c("general", "stephalving"))

  ans           <- .fitmModel(object, method=method, details=details, eps.parm=eps.parm, maxit=maxit, ...)
  ans$dimension <- .mmod_dimension(object$modelinfo$dlq, object$datainfo)

  cg <- object$datainfo$CGstats

  ## Saturated model
  n.obs <- cg$n.obs
  NN    <- cg$N
  qq    <- length(cg$cont.names)
  SS    <- cg$SSD/NN
  logL.sat <- sum(n.obs*log(n.obs/NN)) - NN*qq/2 * log(2*pi) - NN/2 * log(det(SS)) - NN*qq/2

  ## Independence model
  i.model <- loglin(n.obs, as.list(seq_along(dim(n.obs))),iter=1, print=FALSE, fit=TRUE)
  grand.mean <- rowSumsPrim(colwiseProd(n.obs/NN, cg$center))
  SS.ind <- (cg$SS - NN*grand.mean %*% t(grand.mean))/NN
  logL.ind <- sum(n.obs*log(i.model$fit/NN))- NN*qq/2 * log(2*pi) - NN/2 * sum(log(diag(SS.ind))) - NN*qq/2

  ans$ideviance <- ans$logL - logL.ind
  ans$dev       <- logL.sat - ans$logL
  ans$aic       <-  -2*ans$logL + 2*ans$dimension['mod.dim']
  ans$bic       <-  -2*ans$logL + log(nrow(object$datainfo$data))*ans$dimension['mod.dim']

  object$fitinfo  <- ans
  object$isFitted <- TRUE
  return(object)
}

print.MIfit <- function(x,...){
  cat("MIfit:\n")
  cat(sprintf("components: %s \n", toString(names(x))))
  print(x[c("parms","logL","init.logL","dimension")])
  return(invisible(x))
}

#######################################################################
###
### Main fitting function for mModels
###
### 1) Creates initical values
### 2) Finds sufficient statistics
### 3) Calls .mModel_iterate
###
#######################################################################

.fitmModel <- function(object, method="general", details=0, eps.parm=1e-6, eps.logL=1e-6, maxit=100,...){

  t.start <- proc.time()

### Generate initial parameter value
###
  .infoPrint(details,1, "fit.mModel: Creating initial values\n")
  Mparms <- .CGstats2initpms(object$datainfo$CGstats)
  Cparms <- pms2ghkParms(Mparms)
  ##cat("INITIAL PARAMETERS:\n"); print(Mparms)

##   Mp <<- Mparms
##   Cp <<- Cparms
  
### Generator lists
###
  Ad.list <- object$modelinfo$glist.num.disc
  Ac.list <- object$modelinfo$glist.num.cont

### Find weak marginal model for each generator
###
  WMDghk <- .getWeakMarginalDataList(object$datainfo$CGstats,
                                     Ad.list, Ac.list, type="ghk", details=details)
  WMDpms <- lapply(WMDghk, ghk2pmsParms)

##   Wd <<- WMDghk
##   Wp <<- WMDpms
  
### Iterate to maximize logL
###
  .infoPrint(details,1, "fit.mModel: Calling .mModel_iterate\n")
  ans <- .mModel_iterate(Mparms, Cparms, Ad.list, Ac.list,
                         WMDghk, WMDpms, object$datainfo$CGstats,
                         method, eps.parm, eps.logL, maxit, details)

  ans$WMDghk      <- WMDghk
  ans$init.parms  <- Mparms
  class(ans)      <- "MIfit"
  t.use           <- proc.time()-t.start
  #cat("Fitting time:",t.use[1],"\n")
  return(ans)
}

########################################################################
###
### Iteration:
### Repeatedly
### 1) calls .outerloop and
### 2) checks if logL has increased
###
########################################################################

.mModel_iterate <- function(Mparms, Cparms, Ad.list, Ac.list,
                       WMDghk, WMDpms, CGstats,
                       method,
                       eps.parm=1.0e-4, eps.logL=1e-6, maxit=100, details=1){

  prev.Mparms <- Mparms
  init.logL   <- prev.logL <- .mModel_logLpms(CGstats, Mparms)

  .infoPrint(details,3, cat(sprintf("Initial logL %f\n", init.logL)))

  itcount    <- 1
  scale      <- 1

  if (maxit>0){
    repeat {
      zzz       <- .outerloop(Mparms, Cparms, Ad.list, Ac.list,
                              WMDghk, WMDpms, CGstats, scale,
                              method, prev.logL, itcount, details)
      curr.logL <- zzz$logL
      d.logL    <- zzz$d.logL
      d.parms   <- zzz$d.parms

      .infoPrint(details,1,
                 cat(sprintf(".mModel_iterate maxit=%i, curr.logL=%14.6f, d.logL=%12.6f, d.parms=%12.8f\n",
                             maxit, curr.logL, d.logL, d.parms)))

      it.exceed.crit <- itcount >=maxit
      neg.d.logL     <- d.logL < 0

      if (!neg.d.logL & d.logL < eps.logL){
        break() # We are done
      } else {
        if (neg.d.logL){
          cat(sprintf("Fitting method=%s; logL failed to increase (d.logL=%f) - fit may be questionable\n",
                      method,d.logL))
          break()
        } else {
          if (it.exceed.crit){
            cat(sprintf("Fitting method=%s; Maximum number of iterations=%i exceeded - fit may be questionable\n",
                        method,itcount))
            break()
          }
        }
      }

      Mparms      <- zzz$Mparms
      Cparms      <- zzz$Cparms
      prev.Mparms <- Mparms
      prev.logL   <- curr.logL
      itcount     <- itcount + 1
    }
  }
  res <- list(parms=Cparms, logL=curr.logL, init.logL=init.logL)
  res
}

####################################################################
###
### .outerloop
###
### Sets update method; either standard or step-halving
### Iterates over all generators in the model
###
####################################################################

.outerloop <- function(Mparms, Cparms, Ad.list, Ac.list,
                       WMDghk, WMDpms, CGstats, scale,
                       method, logL, itcount, details){

  .infoPrint(details,4, "calling outerloop\n")

  prev.Mparms <- Mparms
  prev.logL   <- logL
  logL.fail   <- 0

  if (method=="general")
    .innerloop <- .standard.innerloop
  else
    .innerloop <- .stephalving.innerloop


  for (ii in 1:length(Ad.list)){
    Ad.idx    <- Ad.list[[ii]]
    Ac.idx    <- Ac.list[[ii]]
    EEghk     <- WMDghk[[ii]]
    EEpms     <- WMDpms[[ii]]
    AApms     <- weakMarginalModel(Mparms, disc=Ad.idx, cont=Ac.idx, type="pms", details=details)
##    AAp <<- AApms
    AAghk     <- pms2ghkParms(AApms)
##    AAg <<- AAghk
    zzz       <- .innerloop(Mparms, Cparms, Ad.idx, Ac.idx,
                            EEghk, EEpms, AAghk, AApms, CGstats, scale, prev.logL, details)
    Mparms    <- zzz$Mparms
    Cparms    <- zzz$Cparms
    curr.logL <- zzz$curr.logL
    logL.fail <- logL.fail + zzz$logL.fail
  }

  curr.logL <- .mModel_logLpms(CGstats,Mparms)
  d.logL    <- curr.logL - prev.logL
  d.parms   <- .mModel_parmdiff(Mparms, prev.Mparms)

  .infoPrint(details,3,
             cat(sprintf("outerloop (%2d): logL %16.10f, d.logL: %16.10f d.parms: %16.10f logL.fail: %f\n",
                         itcount, curr.logL, d.logL, d.parms, logL.fail)))

  list(Mparms=Mparms, Cparms=Cparms, logL=curr.logL, d.logL=d.logL, logL.fail=logL.fail, d.parms=d.parms)
}



.standard.innerloop <- function(Mparms, Cparms, Ad.idx, Ac.idx,
                       EEghk, EEpms, AAghk, AApms, CGstats, scale, prev.logL, details){

  new.Cparms <- .update.ghkParms(Cparms, Ad.idx, Ac.idx,
                                 EEghk, EEpms, AAghk, AApms, scale, CGstats, details=details)

  new.Mparms <- ghk2pmsParms(new.Cparms)
  curr.logL  <- d.logL <- d.parms <- NA
  logL.fail  <- as.numeric(d.logL < 0)

  .infoPrint(details,4,
             cat(sprintf(".std.innerloop(%4.2f): AA=%10s,  curr.logL=%16.10f -2logL=%16.10f d.logL=%16.10f d.parms=%8.6f \n",
                         scale, .toString(c("{",Ad.idx,"|", Ac.idx,"}")), curr.logL, -2*curr.logL, d.logL,
                         d.parms)))

  ans <- list(Mparms=new.Mparms, Cparms=new.Cparms, curr.logL=curr.logL, d.logL=d.logL, logL.fail=logL.fail,
              maxinner.code=NA, step.code=NA)

}


.stephalving.innerloop <- function(Mparms, Cparms, Ad.idx, Ac.idx,
                       EEghk, EEpms, AAghk, AApms, CGstats, scale, prev.logL, details){

  .infoPrint(details,10, "innerloop: finding (model) weak marginals for a generator\n")


  prev.Mparms   <- Mparms
  innercount    <-  1
  maxinner      <-  5
  neg.eps       <-  -1e-4
  good.Mparms   <- Mparms
  good.Cparms   <- Cparms
  step.code     <- 0
  maxinner.code <- 0

  d.logL      <- -99999
  curr.logL   <- prev.logL

  repeat{
    new.Cparms <- .update.ghkParms(good.Cparms, Ad.idx, Ac.idx,
                                   EEghk, EEpms, AAghk, AApms, scale, CGstats, details=details)
    new.Mparms <- ghk2pmsParms(new.Cparms)
    curr.logL  <- .mModel_logLpms(CGstats, new.Mparms)
    d.logL     <- curr.logL - prev.logL
    d.parms    <- .mModel_parmdiff(new.Mparms, prev.Mparms)
    min.eigen  <- min(eigen(new.Mparms$Sigma)$values)

    .infoPrint(details,1,
               cat(sprintf(".steph.innerloop(%4.2f): AA=%10s,  curr.logL=%16.10f -2logL=%16.10f d.logL=%16.10f d.parms=%8.6f \n",
                           scale, .toString(c("{",Ad.idx,"|", Ac.idx,"}")), curr.logL, -2*curr.logL, d.logL,
                           d.parms)))

    if ((d.logL<neg.eps | min.eigen<0) & innercount < maxinner){
      scale       <- scale / 2
      innercount  <- innercount + 1
      step.code   <- 1
    } else {
      if (innercount==maxinner){
        Cparms    <- good.Cparms
        Mparms    <- good.Mparms
        curr.logL <- .mModel_logLpms(CGstats, Mparms)
        maxinner.code <- 1
        .infoPrint(details, 4, cat(sprintf("stephalving failed; restoring original parameters; logL: %10.4f\n", curr.logL)))
      } else {
        Cparms  <- new.Cparms
        Mparms  <- new.Mparms
      }
      break
    }
  }

  logL.fail <- as.numeric(d.logL < 0)
  ans <- list(Mparms=Mparms, Cparms=Cparms, curr.logL=curr.logL, d.logL=d.logL, logL.fail=logL.fail,
              maxinner.code=maxinner.code, step.code=step.code)
  return(ans)
}



.getWeakMarginalDataList <- function(CGstats, Ad.list, Ac.list, type="ghk", details=0){

  .infoPrint(details,1,"fit.mModel: Finding weak (empirical) marginals for each generator\n")

  ##CGstats <- unclass(CGstats)

  ans <- vector("list", length(Ad.list))
  for (ii in 1:length(Ad.list)){
    EE.mm <- weakMarginalData(CGstats, disc=Ad.list[[ii]], cont=Ac.list[[ii]],
                              type=type, details=details)
    ans[[ii]] <- EE.mm
  }
  ans
}

###########################################################################
###
### .update.ghkParms
###
### This is where the parameter updates take place.
###
###########################################################################
.update.ghkParms <- function(Cparms, Ad.idx, Ac.idx, EEghk, EEpms, AAghk, AApms, scale, CGstats, details=0) {

  g.idx <- 1
  h.idx <- 2
  K.idx <- 3

  .infoPrint(details,5, cat(sprintf(".update.ghkParms: A=%8s\n",
                                    .toString(c("{",Ad.idx,"|", Ac.idx,"}")))))
  gt <- .genType(Ad.idx, Ac.idx)
  d.parms.crit <- 0.00001

  if (details>=5){
    cat("PRE UPDATED marginal OBSERVED // FITTED values - moment form\n")
    print(rbind(.as.matrix(ghk2pmsParms(EEghk)),.as.matrix(ghk2pmsParms(AAghk))))
  }

##  cat(".update.ghkParms - calling .mModel_parmdiff\n")
  marg.d.parms <- .mModel_parmdiff(AApms, EEpms)
  .infoPrint(details,5, cat(sprintf("PARMDIF=%f\n", marg.d.parms)))

  if (marg.d.parms>d.parms.crit){
    
    if (details>=5){
      cat("PRE UPDATE Mparms:\n");
      print(.MIparms2matrix(ghk2pmsParms(Cparms)))
      cat("PRE UPDATED marginal OBSERVED // FITTED values - canonical form\n")
      print(rbind(.as.matrix((EEghk)),.as.matrix((AAghk))))
      cat("PRE UPDATE Cparms:\n");
      print(.MIparms2matrix((Cparms)))
    }

    switch(gt,
           "discrete"={
             ##cat("Cparms//Mparms BEFORE update:\n");  print(Cparms)
             upd.g    <- scale*(EEghk[['g']] - AAghk[['g']])
             g.new    <- tableOp2(Cparms[['g']], upd.g, `+`, restore=TRUE)
             #max.chg <- c(max(abs(upd.g)),-1,-1)
             res <- list(g=g.new, h=Cparms[['h']], K=Cparms[['K']])
             res <- .normalize.ghkParms(res)
             res <- c(res[1:3], Cparms[-(1:3)])
             ##cat("Cparms//Mparms AFTER update:\n");  print(res)
           },
           "continuous"={
             h.new   <- Cparms[['h']]
             upd.h   <- scale * (EEghk[['h']]-AAghk[['h']])
             for (jj in 1:ncol(h.new))
               h.new[Ac.idx,jj] <- Cparms[['h']][Ac.idx,jj,drop=FALSE] + upd.h
             upd.k   <- scale * (EEghk[["K"]] - AAghk[["K"]])
             K.new   <- Cparms[['K']]
             K.new[Ac.idx,Ac.idx] <- K.new[Ac.idx,Ac.idx] + upd.k
             ## cat("cont: upd.h:\n"); print(cbind(EEghk[['h']], AAghk[['h']], upd.h))
             ## cat("cont: upd.k:\n"); print(cbind(EEghk[["K"]], AAghk[["K"]], upd.k))
             res <- list(g=Cparms[['g']], h=h.new, K=K.new)
             res <- .normalize.ghkParms(res)
             res <- c(res, Cparms[-(1:3)])
           },
           "mixed"={
             ## g update:
             upd.g    <- scale * (EEghk[[g.idx]] - AAghk[[g.idx]])
             g.new    <- tableOp2(Cparms[[g.idx]], upd.g, `+`, restore=TRUE)

             ##cat("upd.g:\n"); print(t(round(cbind(EEghk[["g"]], AAghk[["g"]], upd.g),4)))
             ## K update:
             upd.k   <- scale * (EEghk[[K.idx]] - AAghk[[K.idx]])
             K.new   <- Cparms[[K.idx]]


             K.new[Ac.idx,Ac.idx] <- K.new[Ac.idx,Ac.idx] + upd.k
             ##cat("upd.k:\n"); print(round(cbind(EEghk[["K"]], AAghk[["K"]], upd.k),4))
             h.new   <- Cparms[[h.idx]]
             upd.h   <- scale * (EEghk[[h.idx]]-AAghk[[h.idx]])

             em       <- AAghk[['jia.mat']]
             Cparms.h <- Cparms[['h']]
             for (jj in 1:ncol(em))
               h.new[Ac.idx,em[,jj]] <- Cparms.h[Ac.idx,em[,jj],drop=FALSE] + upd.h[,jj]
             ##cat("upd.h:\n"); print(round(cbind(EEghk[['h']], AAghk[['h']], upd.h),4))
             ##max.chg <- c(max(abs(upd.g)),max(abs(upd.h)),max(abs(upd.k)))
             res <- list(g=g.new, h=h.new, K=K.new)
             res <- .normalize.ghkParms(res)
             res <- c(res[1:3], Cparms[-(1:3)])
           })
  } else {
    .infoPrint(details, 5, cat(sprintf("Not updating generator\n")))
    res <- Cparms
    #max.chg <- c(0,0,0)
  }

  if (details>=5){
    cat("POST UPDATE Cparms // Mparms:\n");
    MM   <- ghk2pmsParms(Cparms)
    MM$p <- MM$p*MM$N
    RR   <- ghk2pmsParms(res)
    RR$p <- RR$p*RR$N
    print(rbind(.as.matrix(res), .as.matrix(RR)))
  }
  #res$max.chg <- max.chg
  res
}



.mModel_logLpms <- function(CGstats, Mparms){

  Sigma.inv <- solveSPD(Mparms[['Sigma']])

  n.i    <- as.numeric(CGstats[['n.obs']])
  N      <- sum(n.i)
  Q      <- nrow(CGstats[['center']])
  xxx    <- sum(n.i * log(Mparms[['p']])) - N * (Q * log(2*pi) + .logdet(Mparms[['Sigma']])) / 2
  x4     <- -  sum(CGstats[['SSD']] * Sigma.inv) / 2
  mu.dif <- CGstats[['center']] - Mparms[['mu']]
  quad   <- .vMMt(n.i, mu.dif)
  x5     <- - sum(Sigma.inv * quad) / 2
  return(xxx+x4+x5)
}


.mModel_parmdiff <- function(curr.Mparms, prev.Mparms){

##   cat("curr.Mparms:---------------\n "); print(curr.Mparms)
##   cat("prev.Mparms:---------------\n "); print(prev.Mparms)
  
  if (curr.Mparms[['gentype']]=="discrete"){
    N   <- prev.Mparms[['N']]
    cp  <- curr.Mparms[['p']]
    ppp <- as.numeric(N * abs((cp - prev.Mparms[['p']])) /sqrt((N * cp + 1)))
    ans <- max(ppp)
  } else {
    N   <- prev.Mparms[['N']]
    cp  <- curr.Mparms[['p']]
    sss <- curr.Mparms[['Sigma']]
    nr  <- nrow(sss)
    iii <- 1+(nr+1)*((1:nr)-1)
    ddd <- sss[iii] ## faster than diag(sss)
    ##     ppp <- as.numeric(N * abs((cp - prev.Mparms[[p.idx]])) / sqrt((N * cp + 1)))
    ##     mmm <- as.numeric(abs(curr.Mparms[[mu.idx]] - prev.Mparms[[mu.idx]])/sqrt(ddd))

    ppp <- c(N * abs((cp - prev.Mparms[['p']])) / sqrt((N * cp + 1)))
    mmm <- c(abs(curr.Mparms[['mu']] - prev.Mparms[['mu']])/sqrt(ddd))

    xxx <- abs(sss - prev.Mparms[['Sigma']])
    uuu <- xxx/sqrt(tcrossprod(ddd) + sss^2)

    ans <- max(c(ppp,mmm,c(uuu)))
    ##  cat("max parm diff:", ans, "\n")
  }
  ans
}


##############################################################
###
### Create initial parms for mModel from CGstats
###
##############################################################
.CGstats2initpms <- function(CGstats, unif=TRUE){

  ##CGstats <- unclass(CGstats)

  if (unif) { ## Uniform model
    PPP   <- CGstats$n.obs
    PPP[] <- 1 / length(PPP)
    MMM   <- CGstats$center
    MMM[] <- rowMeans(CGstats$center)
    CCC   <- diag(1,nrow(MMM))
  } else {
    ## p.i
    n.i   <- as.numeric(CGstats$n.obs)

    ## mu
    mu.i  <- CGstats$center
    ##mu    <- rowSums(.colmult(n.i, mu.i)) / sum(n.i)
    mu    <- rowSumsPrim(.colmult(n.i, mu.i)) / sum(n.i)

    ## Sigma (total variance when discrete variables are ignored)
    S.i     <- CGstats$cov
    SSD.i   <- .colmult(n.i, S.i)
    ##SSD     <- matrix(rowSums(SSD.i), nrow=nrow(CGstats$center))
    SSD     <- matrix(rowSumsPrim(SSD.i), nrow=nrow(CGstats$center))
    d.mu.i  <- mu.i - mu

    quad    <- .colmult(n.i, d.mu.i) %*% t(d.mu.i)
    #quad <- tcrossprod(.colmult(sqrt(n.i),d.mu.i))
    Sigma   <- (SSD + quad)/sum(n.i)

    ## Create uniform p's
    PPP   <- CGstats$n.obs
    PPP[] <- 1 / length(PPP)

    ## Create uniform means
    MMM   <- CGstats$center
    MMM[] <- rowSumsPrim(MMM)/ncol(MMM)

    CCC   <- Sigma
    if (nrow(CCC)>1)
      CCC <- diag(diag(CCC))
  }

  rownames(CCC) <- colnames(CCC) <- rownames(MMM) <- CGstats$cont.names
  ans        <- c(list(p=PPP,mu=MMM,Sigma=CCC,gentype="mixed"),CGstats[-(1:3)])
  ##class(ans) <- c("pms","MIparms")
  return(ans)
}
