.calc.coeff <- function(Pi, phi, psi, p0, p1, beta0, beta1, dF0, dF1, 
                       t0, t1, tau0, tau1, time.points, interval) {
  calc.time.seq <- function(time.point, times) which(times <= time.point)
  
  q.seq.list <- lapply(time.points, t0, FUN=calc.time.seq)
  r.seq.list <- lapply(time.points, t1, FUN=calc.time.seq)
  names(q.seq.list) <- names(r.seq.list) <- time.points

  tplus0 <- c(pmin(t0, tau0), tau0)
  tplus1 <- c(pmin(t1, tau1), tau1)

  calc.beta.tplus <- function(beta,i,time) return(list(bt=beta*time, i=i))
  beta.tplus0 <- mapply(FUN=calc.beta.tplus,
                        beta=beta0, i=seq_along(beta0),
                        MoreArgs=list(time=tplus0),
                        USE.NAMES=TRUE, SIMPLIFY=FALSE)
  beta.tplus1 <- mapply(FUN=calc.beta.tplus,
                        beta=beta1, i=seq_along(beta1),
                        MoreArgs=list(time=tplus1),
                        USE.NAMES=TRUE, SIMPLIFY=FALSE)

  Pi..p0 <- Pi/p0
  Pi..p1 <- Pi/p1

  coeffs <- mapply(FUN=.calc.Pi.coeff,
                   Pi=Pi, phi=phi, psi=psi,
                   Pi..p0=Pi..p0, Pi..p1=Pi..p1, i=seq.int(along.with=Pi),
                   MoreArgs=list(beta.tplus0=beta.tplus0,
                     beta.tplus1=beta.tplus1,
                     q.seq.list=q.seq.list, r.seq.list=r.seq.list,
                     dF0=dF0, dF1=dF1, interval=interval,
                     tplus0=tplus0, tplus1=tplus1),
                   USE.NAMES = FALSE, SIMPLIFY=FALSE)
  return(coeffs)
}

.calc.Pi.coeff <- function(Pi, phi, psi, tplus0, tplus1,
                           Pi..p0, Pi..p1, i, beta.tplus0, beta.tplus1,
                          q.seq.list, r.seq.list, dF0, dF1, dV, dW, VmF0, WmF1,
                          p0, p1, interval) {
  if(phi == 1)
    return(list(i=i,phi="phi"))
  
  beta0.coeff <- lapply(beta.tplus0, q.list=q.seq.list, dF=dF0, Pi..p=Pi..p0,
                        tplus=tplus0, interval=interval, FUN=.calc.beta.coeff)
  beta1.coeff <- lapply(beta.tplus1, q.list=r.seq.list, dF=dF1, Pi..p=Pi..p1,
                        tplus=tplus1, interval=interval, FUN=.calc.beta.coeff)

  return(list(beta0.coeff=beta0.coeff, beta1.coeff=beta1.coeff, i=i))
}    

.calc.beta.coeff <- function(beta.tplus, q.list, dV, dF, p, Pi..p, Pi.N,
                             tplus, interval) {

  alphahat <- .calc.alphahat(beta.y=beta.tplus$bt, dF=dF, C=Pi..p,
                             interval=interval)

  w <- .calc.w(alpha=alphahat, beta.y=beta.tplus$bt)

  Fas <- sapply(q.list, w.dF=w*dF,
                FUN=function(q.seq, w.dF) sum(w.dF[q.seq])) / Pi..p

  FnAs <- stepfun(x=tplus, c(0, cumsum(w*dF)/Pi..p), right=FALSE)
  return(list(i=beta.tplus$i, alphahat=alphahat, FnAs=FnAs, Fas=Fas))
}

sensitivitySGD <- function(z, s, d, y, v, beta0, beta1, phi, Pi, psi, tau,
                           time.points, selection, trigger, groupings,
                           followup.time,
                           ci=0.95, ci.method=c("bootstrap", "analytic"),
                           ci.type="twoSided",
                           custom.FUN=NULL, na.rm=FALSE, N.boot=100L,
                           N.events=NULL, interval=c(-100,100),
                           upperTest=FALSE, lowerTest=FALSE, twoSidedTest=TRUE,
                           inCore=TRUE, verbose=getOption("verbose"),
                           colsPerFile=1000L, isSlaveMode=FALSE) {

  withoutCdfs <- isSlaveMode && !missing(ci.method) && is.null(ci.method)
  withoutCi <- ((!missing(ci.method) && ci.method == "") ||
                (isSlaveMode && !(!missing(ci.method) &&
                                 !is.null(ci.method) &&
                                 c('analytic') %in% ci.method)))

  doFollowupMethod <- !missing(followup.time) && !is.null(followup.time)
  ## z - group that subject belongs to
  ## s - subject met selection cirteria
  ## d - subject had event
  ## y - time until event ocurred
  if(withoutCdfs)
    ci.method <- NULL
  else if(withoutCi)
    ci.method <- ""
  else
    ci.method <- sort(unique(match.arg(ci.method, several.ok=TRUE)))
  
  
  ErrMsg <- character(0L)
  if(!isSlaveMode) {
    ## Not running a boot strap mode
    ## Run error checks on variables.
    ErrMsg <- NULL
    ErrMsg <- c(.CheckSelection(selection, s),
                .CheckGroupings(groupings),
                .CheckTrigger(trigger, d),
                .CheckTau(tau),
                .CheckPhiPiPsi(phi=phi, Pi=Pi, psi=psi),
                .CheckLength(z=z, s=s, d=d, y=y, v=v),
                .CheckZ(z, groupings, na.rm=na.rm),
                .CheckS(s, na.rm=na.rm),
                .CheckY(y, s, selection, na.rm=na.rm),
                .CheckD(d=d, s=s, selection=selection, na.rm=na.rm),
                .CheckV(v=v, followup.time=followup.time, na.rm=na.rm),
                .CheckCi(ci=ci, ci.type=ci.type))
    
    if(length(ErrMsg) > 0L)
      stop(paste(ErrMsg, collapse="\n  "))

    ## Process tau
    if(length(tau) == 1) {
      tau <- c(tau, tau)
    }

    if(missing(ci.type)) {
      ci.type <- rep('twoSided', length.out=length(ci))
    } else {
      ci.type <- match.arg(ci.type, c('upper', 'lower', 'twoSided'),
                           several.ok=TRUE)
    }

    s <- s == selection

    if(na.rm == TRUE) {
      if(doFollowupMethod) {
        naIndex <- (is.na(s) | is.na(v) | is.na(z) | (s & (is.na(d) | is.na(y))))
      } else {
        naIndex <- (is.na(s) | is.na(z) | (s & (is.na(d) | is.na(y))))
      }

      if(any(naIndex)) {
        z <- z[!naIndex]
        s <- s[!naIndex]
        d <- d[!naIndex]
        y <- y[!naIndex]
        
        if(doFollowupMethod)
          v <- v[!naIndex]
      }
    }

    d <- d == trigger
    
    z <- z == groupings[2L]

  }
 
  test <- c(upper=upperTest, lower=lowerTest, twoSided=twoSidedTest)
  n.test <- sum(test)
  names.test <- names(test)[test]

  ## N  - number subjects
  ## N0 - number of subjects in group 0
  ## N1 - number of subjects in group 1
  N <-length(z)
  N1 <- sum(z)
  N0 <- N-N1

  if(doFollowupMethod) {
    s <- s & v < followup.time

    temp <- with(survfit(Surv(v, s) ~ z, se.fit=FALSE),{
      who <- n.event > 0L
      data.frame(strata=rep.int(seq_along(strata), times=strata)[who],
                 surv=surv[who])
    })

    
    p0 <- 1L - tail(x=temp$surv[temp$strata == 1L], n=1L)
    p1 <- 1L - tail(x=temp$surv[temp$strata == 2L], n=1L)
  } else {
    z0.s1 <- !z & s
    z1.s1 <- z & s
    
    ## n0 - number of subjects in group 0 that were selected 
    ## n1 - number of subjects in group 1 that were selected
    n0 <- sum(z0.s1)
    n1 <- sum(z1.s1)

    ## p0 - probiblity that subject in group 0 was selected
    ## p1 - probablity that subject in group 1 was selected
    p0 <- n0/N0
    p1 <- n1/N1
  }

  if(!isSlaveMode) {
    ErrMsg <- .CheckPhiPiPsi(phi=phi, Pi=Pi, psi=psi, p0=p0, p1=p1)

    if(length(ErrMsg) > 0)
      stop(ErrMsg)
  }
  
  if(all((!z & s) == FALSE) || all((z & s) == FALSE))
    if(isSlaveMode) {
      return(list(SCE = logical(0)))
    } else {
      stop("No events occured in one or more of the treatment arms")
    }

  tmp <- .calcPiPhiPsi(Pi=Pi, phi=phi, psi=psi, p0=p0, p1=p1)
  Pi <- tmp$Pi
  psi <- tmp$psi
  phi <- tmp$phi
  sens.var <- tmp$sens.var
  rm(tmp)

  ## summary survfit of length of time til event for group 0 and group 1.
##   temp <- summary(survfit(Surv(y[z0.s1],d[z0.s1])~1L))
##   t0 <- temp$time
##   F0 <- 1L - temp$surv
   
##   temp <- summary(survfit(Surv(y[z1.s1],d[z1.s1])~1L))
##   t1 <- temp$time
##   F1 <- 1L - temp$surv

  ## Alt implementation
  ## survfit of length of time til event for group 0 and group 1.
  temp <- with(survfit(Surv(y[s], d[s]) ~ z[s], se.fit=FALSE), {
    who <- n.event > 0
    strata <- rep.int(seq_along(strata), times=strata)[who]
    time <- time[who]
    F <- 1L - surv[who]
    data.frame(strata=strata, time=time, F=F)
  })
  
  t0 <- temp$time[temp$strata==1]
  F0 <- temp$F[temp$strata==1]
  t1 <- temp$time[temp$strata==2]
  F1 <- temp$F[temp$strata==2]

  len.t0 <- length(t0)
  len.t1 <- length(t1)

  ## total length of phi vector
  len.total <- 4L + len.t0 + len.t1

  if(len.t0 == 0L || len.t1 == 0) {
    if(isSlaveMode) {
      return(list(SCE = logical(0)))
    } else {
      stop("No times occured in one or more of the treatment arms")
    }
  }

  dF0 <- diff(c(0L,F0,1L))
  dF1 <- diff(c(0L,F1,1L))

  coeffs <- .calc.coeff(Pi=Pi, phi=phi, psi=psi,
                        p0=p0, p1=p1, beta0=beta0, beta1=beta1,
                        dF0=dF0, dF1=dF1, t0=t0, t1=t1, tau0=tau[1],
                        tau1=tau[2], time.points=time.points,
                        interval=interval)

  SCE.dim <- c(length(beta0), length(beta1), length(psi), length(time.points))
  SCE.length <- prod(SCE.dim)
  SCE.dimnames <- list(format(beta0, trim=TRUE),
                       format(beta1, trim=TRUE),
                       format(switch(sens.var,
                                     Pi=Pi,
                                     phi=phi,
                                     psi=psi),
                              trim=TRUE, digits=4,
                              drop0trailing=TRUE),
                       format(time.points, trim=TRUE))
  names(SCE.dimnames) <- c("beta0", "beta1", sens.var, "time.points")

  SCE <- array(numeric(SCE.length), dim=SCE.dim, dimnames=SCE.dimnames)

  if(!is.null(custom.FUN)) {
    result <- array(numeric(SCE.length), dim=SCE.dim, dimnames=SCE.dimnames)
  }

  if(!withoutCdfs) {
    FnAs0.dim <- SCE.dim[c(-2L,-4L)]
    FnAs0.length <- prod(FnAs0.dim)
    FnAs0.dimnames <- SCE.dimnames[c(-2L,-4L)]

    alphahat0 <- array(numeric(FnAs0.length), dim=FnAs0.dim,
                       dimnames=FnAs0.dimnames)

    FnAs1.dim <- SCE.dim[c(-1L, -4L)]
    FnAs1.length <- prod(FnAs1.dim)
    FnAs1.dimnames <- SCE.dimnames[c(-1L, -4L)]

    alphahat1 <- array(numeric(FnAs1.length), dim=FnAs1.dim,
                       dimnames=FnAs1.dimnames)

    FnAs0 <- funArray(vector(mode='list', length=prod(FnAs0.dim)),
                      dim=FnAs0.dim,
                      dimnames=FnAs0.dimnames)

    FnAs1 <- funArray(vector(mode='list', length=prod(FnAs1.dim)),
                      dim=FnAs1.dim,
                      dimnames=FnAs1.dimnames)
  }
  
  ## iterate across all betas and Pi values
  for(Pi.coeff in coeffs) {
    if(is.character(Pi.coeff$phi) && Pi.coeff$phi == "phi") {
      SCE.info <- sensitivitySGL(z=z, s=s, d=d, y=y, beta=beta0,
                                 tau=tau[1L],
                                 time.points=time.points,
                                 selection=selection, trigger=trigger,
                                 groupings=FALSE, ci.method=ci.method,
                                 custom.FUN=custom.FUN, isSlaveMode=TRUE, interval=interval)

      if(!is.null(custom.FUN)) {
        result[,,Pi.coeff$i,] <- SCE.info$result[rep.int(seq_along(beta0),
                                                         times=length(beta1)),]
      }
      if(!withoutCdfs) {
        alphahat0[,Pi.coeff$i] <- SCE.info$alphahat
        FnAs0[,Pi.coeff$i] <- SCE.info$Fas0
        alphahat1[,Pi.coeff$i] <- NA
        FnAs1[,Pi.coeff$i] <- SCE.info$Fas1
      }        
      SCE[,,Pi.coeff$i,] <- SCE.info$SCE[rep.int(seq_along(beta0),
                                                 times=length(beta1)),]
      next
    }
    for(beta0.coeff in Pi.coeff$beta0.coeff) {
      for(beta1.coeff in Pi.coeff$beta1.coeff) {
        SCE[beta0.coeff$i,beta1.coeff$i,Pi.coeff$i,] <- beta0.coeff$Fas - beta1.coeff$Fas
        
        if(!is.null(custom.FUN)) {
          result[beta0.coeff$i,beta1.coeff$i,Pi.coeff$i,] <-
            custom.FUN(Fas0=beta0.coeff$FnAs, Fas1=beta1.coeff$FnAs, time.points=time.points, p0=p0, p1=p1)
        }
      }
      
      if(!withoutCdfs) {
        FnAs0[beta0.coeff$i,Pi.coeff$i] <- beta0.coeff$FnAs
        alphahat0[beta0.coeff$i,Pi.coeff$i] <- beta0.coeff$alphahat
      }
    }

    if(!withoutCdfs) {
      for(beta1.coeff in Pi.coeff$beta1.coeff) {
        FnAs1[beta1.coeff$i,Pi.coeff$i] <- beta1.coeff$FnAs
        alphahat1[beta1.coeff$i,Pi.coeff$i] <- beta1.coeff$alphahat
      }
    }
  }
  
  if(withoutCdfs && !is.null(custom.FUN)) return(list(SCE=SCE, result=result))
  if(withoutCdfs) return(list(SCE = SCE))

  cdfs <- list(alphahat0=alphahat0, beta0=beta0, Fas0=FnAs0,
               alphahat1=alphahat1, beta1=beta1, Fas1=FnAs1,
               psi=psi, phi=phi, Pi=Pi, time.points=time.points)

  if(withoutCi) {
    if(isSlaveMode) {
      if(!is.null(custom.FUN))
        return(c(list(SCE=SCE, result=result), cdfs))
      else
        return(c(list(SCE = SCE), cdfs))
    }

    return(structure(c(list(SCE=SCE,
                            result=if(!is.null(custom.FUN)) result else SCE), cdfs),
                     class=c("sensitivity.1d", "sensitivity"),
                     parameters=list(z0=groupings[1], z1=groupings[2],
                       selected=selection, trigger=trigger)))
  }

  if(!isSlaveMode) {
    ci.map <- vector('list', length(ci.type))
    names(ci.map) <- ci

    for(i in seq_along(ci.type)) {
      if(ci.type[i] == "upper"){
        ci.map[[i]] <- ci[i]
      } else if(ci.type[i] == "lower"){
        ci.map[[i]] <- 1 - ci[i]
      }else if(ci.type[i] == "twoSided") {
        if(ci[i] < 0.5)
          ci.map[[i]] <- c(ci[i] - ci[i]/2, 1 - ci[i]/2)
        else
          ci.map[[i]] <- c((1-ci[i])/2, ci[i] + (1 - ci[i])/2)
      }
    }

    ci.probs <- unique(unlist(ci.map, recursive=FALSE))
    ci.probsLen <- length(ci.probs)

    z.seq <- seq_len(N)

    SCE.ci.dim <- c(SCE.dim, ci.probsLen, length(ci.method))
    SCE.ci.dimnames <- c(SCE.dimnames, list(ci.probs=as.character(ci.probs),
                                            ci.method=ci.method))
    
    SCE.ci <- array(numeric(0),
                    dim=SCE.ci.dim,
                    dimnames=SCE.ci.dimnames)
    if(!is.null(custom.FUN)) {
      result.ci <- array(numeric(0),
                         dim=SCE.ci.dim,
                         dimnames=SCE.ci.dimnames)
    }
  }
  
  SCE.var.dim <- c(SCE.dim, length(ci.method))
  SCE.var.dimnames <- c(SCE.dimnames, list(ci.method=ci.method))
  
  SCE.var <- array(numeric(0),
                   dim=SCE.var.dim,
                   dimnames=SCE.var.dimnames)  

  SCE.p.dim <- c(SCE.dim, n.test, length(ci.method))
  SCE.p.dimnames <- c(SCE.dimnames,
                      list(test=names.test, ci.method=ci.method))

  SCE.p <- array(numeric(0),
                 dim=SCE.p.dim,
                 dimnames=SCE.p.dimnames)
  
  if(!is.null(custom.FUN)) {
    result.var <- array(numeric(0),
                       dim=SCE.var.dim,
                       dimnames=SCE.var.dimnames)
    result.p <- array(numeric(0),
                      dim=SCE.p.dim,
                      dimnames=SCE.p.dimnames)
  }

  if("analytic" %in% ci.method) {
    stop("Analytic method is not currently implemented")
  }

  if(isSlaveMode) {
    if(!is.null(custom.FUN))
      return(c(list(SCE=SCE, SCE.var=SCE.var, result=result), cdfs))
    else
    return(c(list(SCE=SCE, SCE.var=SCE.var), cdfs))
  }
  
  if("bootstrap" %in% ci.method) {
    current.fun <- sys.function()

    N.boot <- as.integer(N.boot)

    if(is.null(N.events)) {
      nVal <- N
      mkBsIndex <- .makeBootstrapLenIndx
    } else {
      nVal <- as.integer(N.events)
      mkBsIndex <- .makeBootstrapEvntIndx
    }

    if(doFollowupMethod) {
      bootCalc <- function(i, z.seq, nVal, beta0, beta1, psi, tau, time.points,
                           current.fun, custom.FUN, interval, verbose) {
        samp <- mkBsIndex(s, indx.seq=z.seq, N=nVal)
        ans <- current.fun(z=z[samp], s=s[samp], v=v[samp],
                           d=d[samp], y=y[samp],
                           beta0=beta0, beta1=beta1, psi=psi,
                           tau=tau, followup.time=followup.time,
                           time.points=time.points, interval=interval,
                           ci.method=NULL, custom.FUN=custom.FUN,
                           isSlaveMode=TRUE)
        if(verbose) cat(".")
        if(!is.null(custom.FUN))
          return(array(c(ans$SCE, ans$result), dim=c(1,length(ans$SCE), 2)))
        else
          return(array(ans$SCE, dim=c(1,length(ans$SCE), 1)))
        return(ans)
      }
    } else {
      bootCalc <- function(i, z.seq, nVal, beta0, beta1, psi, tau, time.points,
                           current.fun, custom.FUN, interval, verbose) {
        samp <- mkBsIndex(s, indx.seq=z.seq, N=nVal)
        ans <- current.fun(z=z[samp], s=s[samp], d=d[samp], y=y[samp],
                           beta0=beta0, beta1=beta1, psi=psi,
                           tau=tau, time.points=time.points,
                           interval=interval, custom.FUN=custom.FUN,
                           ci.method=NULL, isSlaveMode=TRUE)
        
        if(verbose) cat(".")
        if(!is.null(custom.FUN))
          return(array(c(ans$SCE, ans$result), dim=c(1,length(ans$SCE), 2)))
        else
          return(array(ans$SCE, dim=c(1,length(ans$SCE), 1)))
      }
    }
    if(inCore) {
      vals <- do.call(rbind,
                      lapply(integer(N.boot), FUN=bootCalc,
                             z.seq=z.seq, nVal=nVal,
                             beta0=beta0, beta1=beta1,
                             psi=psi, tau=tau,
                             time.points=time.points,
                             current.fun=current.fun,
                             custom.FUN=custom.FUN, interval=interval,
                             verbose=verbose))

      N.bootActual <- nrow(vals)
      
      if(!is.null(custom.FUN))
        dim(vals) <- c(nrow(vals), ncol(vals) %/% 2L, 2L)
      else
        dim(vals) <- c(nrow(vals), ncol(vals), 1L)

      vals <- apply(vals, c(2L,3L),
                    FUN=function(x) return(c(var(x), mean(x > 0), mean(x < 0), quantile(x, probs=ci.probs))))

      SCE.var.boot <- vals[1L,,1L]
      SCE.ci.boot <- t(array(vals[-c(1L,2L,3L),,1L], dim=c(nrow(vals)-3L,ncol(vals))))
      SCE.p.boot <- cbind(if(test['upper']) vals[3,,1L,drop=FALSE],
                          if(test['lower']) vals[2,,1L,drop=FALSE],
                          if(test['twoSided']) 2 * ifelse(vals[2,,1L,drop=FALSE] > vals[3,,1L,drop=FALSE],
                                                          vals[3,,1L,drop=FALSE], vals[2,,1L,drop=FALSE]))

      if(!is.null(custom.FUN)) {
        result.var.boot <- vals[1L,,2L]
        result.ci.boot <- t(array(vals[-c(1L,2L,3L),,2L], dim=c(nrow(vals)-3L, ncol(vals))))
        result.p.boot <- cbind(if(test['upper']) vals[3,,2L,drop=FALSE],
                               if(test['lower']) vals[2,,2L,drop=FALSE],
                               if(test['twoSided']) 2 * ifelse(vals[2,,2L,drop=FALSE] > vals[3,,1L,drop=FALSE],
                                                               vals[3,,2L,drop=FALSE], vals[2,,1L,drop=FALSE]))
      }      
    } else {
      colsPerFile <- as.integer(colsPerFile)

      tmpfile <- tempfile()

      fieldWidth <- nchar(sprintf("%+a", .Machine$double.xmax))
      recordWidth <- fieldWidth + 1L
      lineWidth <- recordWidth*SCE.length
      fieldFmt <- sprintf("%%+%da%s", fieldWidth,
                          rep(c(" ", "\n"), times=c(SCE.length - 1L, 1L)))
      
      lapply(logical(N.boot),
             FUN=function(...) cat(sprintf(fieldFmt, bootCalc(...)), sep="",
               file=tmpfile, append=TRUE),
             z.seq=z.seq, nVal=nVal,
             beta0=beta0, beta1=beta1,
             psi=psi, tau=tau,
             time.points=time.points,
             current.fun=current.fun,
             verbose=verbose)
      if(verbose) cat('\n')

      
      needFiles <- (SCE.length %/% colsPerFile)
      remainder <- SCE.length %% colsPerFile

      filesNCols <- rep(c(colsPerFile, remainder), times=c(needFiles, remainder > 0L))
      readWidths <- filesNCols*recordWidth
      outFilenames <- sprintf("%s_split_%0*d", tmpfile,
                              nchar(length(readWidths)),
                              seq_along(readWidths))
      inConn <- file(tmpfile, open="r")
      on.exit(close(inConn))

      cullColumns <- function(offset, readWidth, outFilename, inConn,
                              lineWidth, N.boot) {
        outConn <- file(outFilename, open='w')
        on.exit(close(outConn))
        seek(inConn, where=offset, origin="start")

        lapply(logical(N.boot),
               FUN = function(j, inConn, outConn, readWidth, lineWidth) {
                 section <- readChar(inConn, readWidth)

                 if(substr(section, readWidth, readWidth) == " ")
                   substr(section, readWidth, readWidth) <- "\n"

                 writeChar(section, outConn, readWidth, eos=NULL)

                 seek(inConn, where=lineWidth - readWidth, origin="current")
                 return(NULL)
               }, inConn=inConn, outConn=outConn, readWidth=readWidth, lineWidth=lineWidth)

        cat("*")
        return(NULL)
      }
      
      mapply(FUN=cullColumns,
             offset=cumsum(c(0L, readWidths[-1L])),
             readWidth=readWidths,
             outFilename=outFilenames,
             MoreArgs=list(inConn=inConn, lineWidth=lineWidth, N.boot=N.boot))
      if(verbose) cat("\n")
      
      SCE.boot <- do.call(cbind, mapply(MoreArgs=list(ci.probs=ci.probs),
                                        FUN=function(filename, numCols, ci.probs) {
                                          dat <- scan(filename,
                                                      what=rep(list(double(0L)), numCols),
                                                      quiet=TRUE)
                                          ans <- sapply(dat,
                                                        FUN=function(column) c(var(column), quantile(column, ci.probs, names=TRUE)))
                                          
                                          if(verbose) cat('.')

                                          return(ans)
                                        },
                                        filename=outFilenames, numCols=filesNCols))
      if(verbose) cat('\n')
    }

    dim(SCE.var.boot) <- SCE.dim
    SCE.var[,,,,"bootstrap"] <- SCE.var.boot

    dim(SCE.ci.boot) <- c(SCE.dim, ci.probsLen)
    SCE.ci[,,,,,"bootstrap"] <- SCE.ci.boot

    SCE.p[,,,,,"bootstrap"] <- SCE.p.boot

    if(!is.null(custom.FUN)) {
      dim(result.var.boot) <- SCE.dim
      result.var[,,,,"bootstrap"] <- result.var.boot

      dim(result.ci.boot) <- c(SCE.dim, ci.probsLen)
      result.ci[,,,,,"bootstrap"] <- result.ci.boot

      SCE.p[,,,,,"bootstrap"] <- result.p.boot
    }
  }

  ans <- structure(c(list(SCE=SCE, SCE.ci=SCE.ci, SCE.var=SCE.var, SCE.p=SCE.p),
                     if(!is.null(custom.FUN))
                       list(result=result, result.ci=result.ci,
                            result.var=result.var, result.p=result.p),
                     cdfs, list(ci.map=ci.map)),
                   class=c("sensitivity.1d", "sensitivity"),
                   parameters=list(z0=groupings[1], z1=groupings[2],
                     selected=selection, trigger=trigger))

  if('bootstrap' %in% ci.method) {
    attr(ans, 'N.boot') <- N.boot
    attr(ans, 'N.bootActual') <- N.bootActual
  }

  return(ans)
}
