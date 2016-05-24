arx <-
function(y, mc=FALSE, ar=NULL, ewma=NULL, mxreg=NULL,
  vc=FALSE, arch=NULL, asym=NULL, log.ewma=NULL, vxreg=NULL,
  zero.adj=0.1, vc.adj=TRUE,
  vcov.type=c("ordinary", "white", "newey-west"),
  qstat.options=NULL, tol=1e-07, LAPACK=FALSE, verbose=TRUE,
  plot=TRUE)
{
  vcov.type <- match.arg(vcov.type)

  ##regressand:
  y.name <- deparse(substitute(y))
  if(is.zoo(y)){ y <- cbind(y) }else{ y <- as.zoo(cbind(y)) }
  #OLD: y <- as.zoo(cbind(y))
  y <- cbind(y)
  if(NCOL(y) > 1) stop("Dependent variable not 1-dimensional")
  if( is.null(y.name)){ y.name <- colnames(y)[1] }
  if( y.name[1] =="" ){ y.name <- "y" }
  y.n <- NROW(y)
  y.index <- index(y)
  y <- coredata(y)

#OLD:
#  y <- na.trim(y)
#  y.n <- NROW(y)
#  y.index <- index(y)
#  t1 <- y.index[1]
#  t2 <- y.index[y.n]
#  y <- coredata(y)
#  y <- y[,1]

  ##regressors:
  mX <- NULL
  mXnames <- NULL

  ##mean intercept:
  if(as.numeric(mc)==1){
    mX <- cbind(rep(1,y.n))
    mXnames  <- "mconst"
  }

  ##ar terms:
  if(!is.null(ar)){
    tmp <- NULL
    nas <- rep(NA, max(ar))
    tmpfun <- function(i){
      tmp <<- cbind(tmp, c(nas[1:i],y[1:c(y.n-i)]))
    }
    tmpfun <- sapply(ar,tmpfun)
    mX <- cbind(mX, tmp)
    mXnames <- c(mXnames, paste("ar", ar, sep=""))
  }

  ##ewma term:
  if(!is.null(ewma)){
    tmp <- do.call(eqwma, c(list(y),ewma) )
    mXnames <- c(mXnames, colnames(tmp))
    colnames(tmp) <- NULL
    mX <- cbind(mX, tmp)
  }

  ##adjust for NAs:
  tmp <- zoo(cbind(y,mX), order.by=y.index)
  tmp <- na.trim(tmp, sides="both", is.na="any")
  y <- tmp[,1]
  y.n <- NROW(y) #re-define y.n
  y.index <- index(y) #re-define y.index
  t1 <- y.index[1]
  t2 <- y.index[y.n]
  y <- coredata(y)
  if(!is.null(mX)){
    mX <- tmp[,2:NCOL(tmp)]
    mX <- coredata(mX)
    mX <- cbind(mX)
    colnames(mX) <- NULL
  }

  ##mxreg:
  if(!is.null(mxreg)){
    mxreg <- as.zoo(cbind(mxreg))
    mxreg.names <- colnames(mxreg)
    if(is.null(mxreg.names)){
      mxreg.names <- paste("mxreg", 1:NCOL(mxreg), sep="")
    }
    if(any(mxreg.names == "")){
      missing.colnames <- which(mxreg.names == "")
      for(i in 1:length(missing.colnames)){
        mxreg.names[i] <- paste("mxreg", i, sep="")
      }
    }
    mXnames <- c(mXnames, mxreg.names)
    mxreg <- window(mxreg, start=t1, end=t2)
    mxreg <- cbind(coredata(mxreg))
    mX <- cbind(mX, mxreg)

    ##re-adjust for NAs:
    tmp <- zoo(cbind(y,mX), order.by=y.index)
    tmp <- na.trim(tmp, sides="both", is.na="any")
    y <- tmp[,1]
    y.n <- NROW(y) #re-define y.n
    y.index <- index(y) #re-define y.index
    t1 <- y.index[1]
    t2 <- y.index[y.n]
    y <- coredata(y)
    mX <- tmp[,2:NCOL(tmp)]
    mX <- coredata(mX)
    mX <- cbind(mX)
    colnames(mX) <- NULL

  } #end if(!is.null(mxreg))

  ##vxreg:
  if(!is.null(vxreg)){
    vxreg <- as.zoo(cbind(vxreg))
    vxreg.names <- colnames(vxreg)
    if(is.null(vxreg.names)){
      vxreg.names <- paste("vxreg", 1:NCOL(vxreg), sep="")
    }
    if(any(vxreg.names == "")){
      missing.colnames <- which(vxreg.names == "")
      for(i in 1:length(missing.colnames)){
        vxreg.names[i] <- paste("vxreg", i, sep="")
      }
    }
    vxreg <- window(vxreg, start=t1, end=t2)
#OLD:
#    vxreg <- cbind(coredata(vxreg))
    colnames(vxreg) <- NULL
  } #end if(!is.null(vxreg))

  ##determine qstat.options:
  if(is.null(qstat.options)){
    if(is.null(ar)){ar.lag <- 1}else{ar.lag <- max(ar)+1}
    if(is.null(arch)){arch.lag <- 1}else{arch.lag <- max(arch)+1}
    qstat.options <- c(ar.lag, arch.lag)
  }

#OLD:
#  ##adjust y, regressors and index:
#  ewma.mean.chk <- if(is.null(ewma)){0}else{ifelse(is.null(ewma$lag),1,ewma$lag)}
#  max.ar <- if( is.null(ar)&&is.null(ewma) ){0}else{ max(ar,ewma.mean.chk) }
#  if(max.ar > 0){
#    y <- y[c(max.ar+1):y.n]
#    y.index <- y.index[c(max.ar+1):y.n]
#    if(!is.null(mX)){mX <- cbind(mX[c(max.ar+1):y.n,])}
#    if(!is.null(vxreg)){
#      vxreg <- cbind(vxreg[c(max.ar+1):y.n,])
#    }
#    y.n <- length(y) #new length
#  }

  ##aux: info for getsm/getsv functions
  aux <- list()
  aux$y <- y
  aux$y.index <- y.index
  aux$y.name <- y.name
  aux$y.n <- y.n
  if(!is.null(mX)){
    colnames(mX) <- NULL
    aux$mX <- mX
    aux$mXnames <- mXnames
    aux$mXncol <- NCOL(mX)
  }
  aux$vc <- vc
  aux$zero.adj <- zero.adj
  aux$vc.adj <- vc.adj
  aux$vcov.type <- vcov.type
  aux$qstat.options <- qstat.options
  aux$tol <- tol
  aux$LAPACK <- LAPACK

  ### INITIALISE ##########

  out <- list()
  out$call <- sys.call()

  #### MEAN ###############

  if(is.null(mX)){

    resids <- aux$y
    fit.m <- rep(0,aux$y.n)
    mean.results <- NULL
    vcov.mean <- NULL

  }else{

    ##estimates, fitted, residuals, etc.:
    est.m <- ols(y, mX, tol = tol, LAPACK=LAPACK,
      method=2)
    fit.m <- as.vector(mX%*%cbind(est.m$coefficients))
    resids <- y - fit.m
    resids2 <- resids^2
    d.f. <- y.n - aux$mXncol
    sigma2 <- sum(resids2)/d.f.

    ##estimate s.e.; compute t-stats. and p-vals.:
    if(vcov.type == "ordinary"){
      varcovmat <- sigma2*est.m$xtxinv
      s.e. <- sqrt(as.vector(diag(varcovmat)))
    }
    if(vcov.type == "white"){
      omega.hat <- crossprod(mX, mX*resids2)
      varcovmat <- est.m$xtxinv %*% omega.hat %*% est.m$xtxinv
      coef.var <- as.vector(diag(varcovmat))
      s.e. <- sqrt(coef.var)
    }
    if(vcov.type == "newey-west"){
      iL <- round(y.n^(1/4), digits=0)
      vW <- 1 - 1:iL/(iL+1)
      vWsqrt <- sqrt(vW)
      mXadj <- resids*mX
      mS0 <- crossprod(mXadj)

      mSum <- 0
      for(l in 1:iL){
        mXadjw <- mXadj*vWsqrt[l]
        mXadjwNo1 <- mXadjw[-c(1:l),]
        mXadjwNo2 <- mXadjw[-c(c(y.n-l+1):y.n),]
        mSum <- mSum + crossprod(mXadjwNo1, mXadjwNo2) + crossprod(mXadjwNo2, mXadjwNo1)
      }

      omega.hat <- mS0 + mSum
      varcovmat <- est.m$xtxinv %*% omega.hat %*% est.m$xtxinv
      coef.var <- as.vector(diag(varcovmat))
      s.e. <- sqrt(coef.var)
    }
    colnames(varcovmat) <- mXnames
    rownames(varcovmat) <- mXnames
    vcov.mean <- varcovmat
    t.stat <- est.m$coefficients/s.e.
    p.val <- pt(abs(t.stat), d.f., lower.tail=FALSE)*2

    mean.results <- as.data.frame(cbind(est.m$coefficients, s.e., t.stat, p.val))
    colnames(mean.results) <- c("coef", "std.error", "t-stat", "p-value")
    rownames(mean.results) <- mXnames

  } #end if(is.null(mX))else(..)

  #### VARIANCE #############

  ##check if log-arch spec:
  logvar.spec.chk <- if( vc==FALSE && is.null(arch)
    && is.null(asym) && is.null(log.ewma)
    && is.null(vxreg) ){0}else{1}

  ##if no log-arch spec:
  if( logvar.spec.chk==0 ){
    loge2.n <- y.n
    aux$loge2.n <- loge2.n
    fit.v <- var(resids)
    resids.std <- resids/sqrt(fit.v)
    vcov.var <- NULL
    variance.results <- NULL
  }

  ##if log-arch spec:
  if( logvar.spec.chk==1 ){

    ##regressand
    zero.where <- which(resids==0)
    eabs <- abs(resids)
    if(length(zero.where) > 0){
      eabs[zero.where] <- quantile(eabs[-zero.where], zero.adj)
    }
    loge2 <- log(eabs^2)

    ##regressor matrix:
    vX <- cbind(rep(1,y.n))
    vXnames <- "vconst"

    ##arch terms:
    if(!is.null(arch)){
      tmp <- NULL
      nas <- rep(NA, max(arch))
      tmpfun <- function(i){
        tmp <<- cbind(tmp, c(nas[1:i],loge2[1:c(y.n-i)]))
      }
      tmpfun <- sapply(arch,tmpfun)
      vX <- cbind(vX, tmp)
      vXnames <- c(vXnames, paste("arch", arch, sep=""))
    }

    ##asym terms:
    if(!is.null(asym)){
      tmp <- NULL
      nas <- rep(NA, max(asym))
      tmpfun <- function(i){
        tmp <<- cbind(tmp, c(nas[1:i],
          loge2[1:c(y.n-i)]*as.numeric(resids[1:c(y.n-i)]<0)))
      }
      tmpfun <- sapply(asym,tmpfun)
      vX <- cbind(vX, tmp)
      vXnames <- c(vXnames, paste("asym", asym, sep=""))
    }

    ##log.ewma term:
    if(!is.null(log.ewma)){
      if(is.list(log.ewma)){
        log.ewma$lag <- 1
      }else{
        log.ewma <- list(length=log.ewma)
      }
      tmp <- do.call(leqwma, c(list(resids),log.ewma) )
      vXnames <- c(vXnames, colnames(tmp))
      colnames(tmp) <- NULL
      vX <- cbind(vX, tmp)
    }

    ##adjust for NAs:
    tmp <- zoo(cbind(loge2,vX), order.by=y.index)
    tmp <- na.trim(tmp, sides="left", is.na="any")
    loge2 <- tmp[,1]
    loge2.n <- NROW(loge2)
    loge2.index <- index(loge2) #re-define y.index
    loge2 <- coredata(loge2)
    vX <- tmp[,2:NCOL(tmp)]
    vX <- coredata(vX)
    vX <- cbind(vX)
    colnames(vX) <- NULL

    ##vxreg:
    if(!is.null(vxreg)){
      vxreg <- window(vxreg, start=loge2.index[1],
        end=loge2.index[loge2.n])
      vxreg <- cbind(vxreg)
      vX <- cbind(vX, coredata(vxreg))
      vXnames <- c(vXnames, vxreg.names)
      colnames(vxreg) <- vxreg.names

      ##re-adjust for NAs:
      tmp <- zoo(cbind(loge2,vX), order.by=loge2.index)
      tmp <- na.trim(tmp, sides="left", is.na="any")
      loge2 <- tmp[,1]
      loge2.n <- NROW(loge2)
      loge2.index <- index(loge2) #re-define index
      loge2 <- coredata(loge2)
      vX <- tmp[,2:NCOL(tmp)]
      vX <- coredata(vX)
      vX <- cbind(vX)
      colnames(vX) <- NULL
    }

#OLD:
#    ##adjust loge2 and regressors:
#    ewma.vol.chk <- if(is.null(log.ewma)){0}else{
#     ifelse(is.null(log.ewma$lag), 1, log.ewma$lag)
#    }
#    pstar <- max(arch, asym, ewma.vol.chk)
#    if(pstar > 0){
#      loge2 <- loge2[c(pstar+1):y.n]
#      vX <- cbind(vX[c(pstar+1):y.n,])
#      vX.index <- y.index[c(pstar+1):y.n]
#    }
#    loge2.n <- length(loge2)

    ##aux: more info for getsm/getsv functions
    aux$loge2 <- loge2
    aux$loge2.n <- loge2.n
    aux$vX <- vX
    aux$vXnames <- vXnames
    aux$vXncol <- NCOL(vX)
    aux$arch <- arch
    aux$asym <- asym
    aux$log.ewma <- log.ewma
    aux$vxreg <- vxreg #note: NROW(vxreg)!=NROW(vX) is possible

    ##estimate:
    est.v <- ols(loge2, vX, tol=tol, LAPACK=LAPACK,
      method=2)
    fit.v <- as.vector(vX%*%cbind(est.v$coefficients))
    ustar <- loge2-fit.v
    ustar2 <- ustar^2
    d.f.v. <- loge2.n - aux$vXncol
    sigma2.v <- sum(ustar2)/d.f.v.

    ##covariance coefficient matrix:
    varcovmat.v <- sigma2.v*est.v$xtxinv
    s.e. <- sqrt(as.vector(diag(varcovmat.v)))
    colnames(varcovmat.v) <- vXnames
    rownames(varcovmat.v) <- vXnames
    vcov.var <- varcovmat.v[-1,-1]
    t.stat <- est.v$coefficients/s.e.
    p.val <- pt(abs(t.stat), d.f.v., lower.tail=FALSE)*2

    if(vc.adj){
      Elnz2 <- -log(mean(exp(ustar)))
      t.stat[1] <- ((est.v$coefficients[1]-Elnz2)^2)/s.e.[1]^2
      p.val[1] <- pchisq(t.stat[1], 1, lower.tail=FALSE)
      est.v$coefficients[1] <- est.v$coefficients[1] - Elnz2
    }
    fit.v <- exp(fit.v - Elnz2)
    resids.std <- resids[c(y.n-loge2.n+1):y.n]/sqrt(fit.v)
#OLD:
#    resids.std <- resids[c(pstar+1):y.n]/sqrt(fit.v)

    variance.results <- as.data.frame(cbind(est.v$coefficients, s.e., t.stat, p.val))
    colnames(variance.results) <- c("coef", "std.error", "t-stat", "p-value")
    rownames(variance.results) <- vXnames
  }

  ### DIAGNOSTICS #################

  if(verbose){
    diagnostics.table <- matrix(NA, 4, 3)
    colnames(diagnostics.table) <- c("Chi-sq", "df", "p-value")
    rownames(diagnostics.table) <- c(paste("Ljung-Box AR(", qstat.options[1], ")", sep=""),
      paste("Ljung-Box ARCH(", qstat.options[2], ")", sep=""),
      "Jarque-Bera", "R-squared")
    ar.LjungBox <- Box.test(resids.std, lag = qstat.options[1], type="L")
    diagnostics.table[1,1] <- ar.LjungBox$statistic
    diagnostics.table[1,2] <- qstat.options[1]
    diagnostics.table[1,3] <- ar.LjungBox$p.value
    arch.LjungBox <- Box.test(resids.std^2, lag = qstat.options[2], type="L")
    diagnostics.table[2,1] <- arch.LjungBox$statistic
    diagnostics.table[2,2] <- qstat.options[2]
    diagnostics.table[2,3] <- arch.LjungBox$p.value
    normality.test <- jb.test(resids.std)
    diagnostics.table[3,1] <- normality.test$statistic
    diagnostics.table[3,2] <- 2
    diagnostics.table[3,3] <- normality.test$p.value

    ##R-squared:
    TSS <- sum( (y - mean(y))^2 )
    RSS <- sum( (resids - mean(resids))^2 )
    Rsquared <- 1 - RSS/TSS
    diagnostics.table[4,1] <- Rsquared
    out$diagnostics <- diagnostics.table
  } #end if(verbose)

  ### OUTPUT: ######################

  out$resids <- zoo(resids, order.by=y.index)
  if(verbose){
    ##mean
    out$mean.fit <- zoo(fit.m, order.by=y.index)

    ##log-variance
    add.nas2var <- rep(NA,y.n-loge2.n)
    out$var.fit <- zoo(c(add.nas2var,fit.v),
      order.by=y.index)
    if(logvar.spec.chk==1){
      out$resids.ustar <- zoo(c(add.nas2var,ustar),
        order.by=y.index)
    }
    out$resids.std <- zoo(c(add.nas2var,resids.std),
      order.by=y.index)
    if(vc.adj && logvar.spec.chk==1){ out$Elnz2 <- Elnz2 }
#DELETE?:
    out$logl <- -loge2.n*log(2*pi)/2 - sum(log(fit.v))/2 - sum(resids.std^2)/2
  } #end if(verbose)

  ##result:
  out$vcov.mean <- vcov.mean
  out$vcov.var <- vcov.var
  out$mean.results <- mean.results
  out$variance.results <- variance.results
  out <- c(list(date=date(),aux=aux), out)
  class(out) <- "arx"
  if(plot){ plot.arx(out) }
  return(out)
}
