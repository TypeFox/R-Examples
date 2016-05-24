sensitivityGBH <- function(z, s, y, beta, selection, groupings,
                           empty.principal.stratum, ci=0.95,
                           ci.method=c("analytic", "bootstrap"),
                           ci.type="twoSided", custom.FUN=NULL, na.rm=FALSE, 
                           N.boot=100, interval=c(-100,100),
                           upperTest=FALSE, lowerTest=FALSE, twoSidedTest=TRUE,
                           method=c("ACE", "T1", "T2"),
                           isSlaveMode=FALSE)
{
  withoutCdfs <- isSlaveMode && !missing(ci.method) && is.null(ci.method)
  withoutCi <- ((!isSlaveMode && !missing(ci.method) && ci.method == "") ||
                (isSlaveMode && !(!missing(ci.method) &&
                                 !is.null(ci.method) &&
                                 'analytic' %in% ci.method)))

  doInfinite <- any(is.infinite(beta))
  doFinite <- any(is.finite(beta))
  hasCustomFun <- !is.null(custom.FUN)
  

  ## if(doInfinite && !doFinite) {
  ##   funCall <- match.call()
  ##   funCall$beta <- c("upper", "lower")[match(beta, c(-Inf, Inf), nomatch=0)]
  ##   names(funCall)[names(funCall) == "beta"] <- "bound"
  ##   funCall[[1]] <- as.name("sensitivityHHS")

  ##   return(eval(funCall, envir=parent.frame()))
  ## }

  if(!isSlaveMode) {
    ## Not running a boot strap mode
    ## Run error checks on variables.
    ErrMsg <- c(.CheckEmptyPrincipalStratum(empty.principal.stratum),
                .CheckSelection(selection, s, empty.principal.stratum),
                .CheckGroupings(groupings),
                .CheckLength(z=z, s=s, y=y),
                .CheckZ(z, groupings, na.rm=na.rm),
                .CheckS(s, empty.principal.stratum, na.rm=na.rm),
                .CheckY(y, s, selection, na.rm=na.rm),
                .CheckCi(ci=ci, ci.type=ci.type))

    if(length(ErrMsg) > 0L)
      stop(paste(ErrMsg, collapse="\n  "))
    rm(ErrMsg)

    s <- s == selection

    if(na.rm == TRUE) {
      naIndex <- !(is.na(s) | is.na(z) | (s & is.na(y)))

      z <- z[naIndex]
      s <- s[naIndex]
      y <- y[naIndex]
    }
    
    if(missing(ci.type)) {
      ci.type <- rep('twoSided', length.out=length(ci))
    } else {
      ci.type <- match.arg(ci.type, c('upper', 'lower', 'twoSided'),
                           several.ok=TRUE)
    }

    GroupReverse <- FALSE
    if(empty.principal.stratum[1] == selection) {
      z <- z == groupings[1]
      GroupReverse <- TRUE
    } else if(empty.principal.stratum[2] == selection)
      z <- z == groupings[2]
  } else {
    GroupReverse <- groupings
  }

  method <- match.bitarg(method)
  n.method <- sum(method)
  names.method <- names(method)[method]
  
  test <- c(upper=upperTest, lower=lowerTest, twoSided=twoSidedTest)
  n.test <- sum(test)
  names.test <- names(test)[test]
  
  if(withoutCi)
    ci.method <- NULL
  else if(isSlaveMode)
    ci.method <- "analytic"
  else
    ci.method <- sort(unique(match.arg(ci.method, several.ok=TRUE)))

  n.ci.method <- length(ci.method)

  if(any(is.na(z) | is.na(s)))
    stop("s, z cannot contain any NA values")
  
  if(any(s & is.na(y)))
    stop("selected y values cannot be NA")
  
  beta.orig <- beta
  beta <- unique(beta)
  bIndex <- match(beta.orig, beta)

  finiteBeta <- beta[is.finite(beta)]
  finiteIndex <- match(finiteBeta, beta, nomatch=0)
  infiniteBeta <- beta[is.infinite(beta)]
  infiniteIndex <- match(infiniteBeta, beta, nomatch=0)

  bounds <- c("upper", "lower")[match(infiniteBeta, c(-Inf, Inf), nomatch=0)]

  z0.s1 <- !z & s
  z1.s1 <- z & s
  z.s1 <- z[s]
  
  y0 <- y[z0.s1]
  y1 <- y[z1.s1]
  
  N <- length(z)
  N0 <- sum(!z)
  n0 <- sum(z0.s1)
  p0 <- n0/N0

  N1 <- sum(z)
  n1 <- sum(z1.s1)
  p1 <- n1/N1

  RR <- min(p1/p0, 1)
  
  R <- rank(y[s])

  Ras0 <- ecdf(R[!z.s1])
  R0.uniq <- knots(Ras0)
  dR0 <- diff(c(0, Ras0(R0.uniq)))
  
  Fas0 <- ecdf(y0)
  y0.uniq <- knots(Fas0)
  dF0 <- diff(c(0, Fas0(y0.uniq)))
  
  Fas1 <- ecdf(y[z1.s1])
  y1.uniq <- knots(Fas1)
  dFas1 <- diff(c(0, Fas1(y1.uniq)))

  y0.mean <- mean(y[z0.s1], na.rm=TRUE)
  y1.mean <- mean(y[z1.s1], na.rm=TRUE)
  
  Rmu1 <- mean(R[z.s1])

  calc.alphaAndW <- function(beta, y, dF, C, interval) {
    alphahat <- .calc.alphahat(beta.y=beta*y, dF=dF, C=C, interval=interval)
    
    w <- .calc.w(alpha=alphahat, beta.y=beta*y)
    return(c(alphahat, w))
  }

  calc.mu0AndAlpha <- function(betas, y, R, dF, dR, C, y0, y1, Rmu1,
                                  n0, n1, y0.mean, interval, method, GroupReverse) {
    method["ACE"] <- TRUE
    resp <- vector(mode="list", length=sum(method) + 2L + !is.null(custom.FUN))
    names(resp) <- c("alphahat", "dFas", names(method)[method],
                     if(!is.null(custom.FUN)) "result")

    temp <- matrix(sapply(betas, FUN=calc.alphaAndW, y=y, dF=dF, C=C,
                          interval=interval), ncol=length(betas))

    resp$alphahat <- temp[1L,]
    w <- temp[-1L,, drop=FALSE]

    # Done with temp
    rm(temp)
    
    resp$dFas <- dF*w/C

    resp$mu0 <- colSums(y*resp$dFas)
    if(method["T1"]) {
      Rmu0 <- colSums(R*(dR0*w/C))
      if(GroupReverse)
        resp$T1 <- Rmu0 - Rmu1
      else
        resp$T1 <- Rmu1 - Rmu0

      ## Done with Rmu0
    }

    if(method["T2"]) {
      Y0star <- outer(y, y0.mean - resp$mu0,
                      FUN=`-`)[rep.int(seq_along(y),
                                       times=tabulate(match(y0, y))),,
                               drop=FALSE]
      
      resp$T2 <- sapply(split(Y0star, col(Y0star)),
                        FUN=function(Y0star, Y1star, z) {
                          Rstar <- rank(c(Y0star, Y1star))

                          if (GroupReverse)
                            return(mean(Rstar[!z]) - mean(Rstar[z]))
                          else
                            return(mean(Rstar[z]) - mean(Rstar[!z]))
                        }, Y1star=y1,
                        z=rep.int(c(FALSE,TRUE), times=c(n0,n1)))

      ## Done with Y0star
    }

    return(resp)
  }

  if(!withoutCdfs)
    Fas0 <- funVector(length(beta))
  
  alphahat <- numeric(length(beta))
  ACE.dim <- length(beta)
  ACE.length <- prod(ACE.dim)
  ACE.dimnames <- format(beta, trim=TRUE)

  temp <- numeric(ACE.dim)
  names(temp) <- ACE.dimnames

  if(method["ACE"])
    ACE <- temp

  if(method["T1"])
    T1 <- temp

  if(method["T2"])
    T2 <- temp

  if(hasCustomFun)
    result <- temp
  
  #Done with temp
  rm(temp)
  
  if(doFinite) {
    temp <- calc.mu0AndAlpha(finiteBeta, y=y0.uniq, R=R0.uniq, dF=dF0, dR=dR0,
                             C=RR, y0=y0, y1=y1, Rmu1=Rmu1, n0=n0, n1=n1,
                             y0.mean=y0.mean, interval=interval, method=method,
                             GroupReverse=GroupReverse)
    
    alphahat[finiteIndex] <- temp$alphahat
    dFas0 <- temp$dFas

    mu0 <- temp$mu0

    if(method["T1"])
      T1[finiteIndex] <- temp$T1

    if(method["T2"])
      T2[finiteIndex] <- temp$T2

    ## Done with temp
    rm(temp)
    
    if(!withoutCdfs) {
      Fas0[finiteIndex] <- lapply(X=as.data.frame(dFas0), FUN=function(x) {
        x <- cumsum(x)
        stepfun(y0.uniq, y=c(0, x))
      })
    }

    if(method["ACE"]) {
      if(GroupReverse) {
        ACE[finiteIndex] <- mu0 - y1.mean
      } else {
        ACE[finiteIndex] <- y1.mean - mu0
      }
    }

    if(hasCustomFun)
      result[finiteIndex] <- custom.FUN(mu0=mu0, mu1=y1.mean, p0=p0, p1=p1)
  }

  if(doInfinite) {
    if(withoutCdfs) {
      infiniteACE.info <- sensitivityHHS(z=z,s=s,y=y, bound=bounds,
                                         groupings=GroupReverse,
                                         ci.method=NULL, custom.FUN=custom.FUN,
                                         method=names.method, isSlaveMode=TRUE)
    } else {
      infiniteACE.info <- sensitivityHHS(z=z,s=s,y=y, bound=bounds,
                                         groupings=GroupReverse,
                                         ci.method=ci.method,
                                         custom.FUN=custom.FUN,
                                         method=names.method, isSlaveMode=TRUE)

      if(!withoutCdfs) {
        Fas0[infiniteIndex] <- infiniteACE.info$Fas0
      }
      
      alphahat[infiniteIndex] <- NA
    }

    if(method["ACE"])
      ACE[infiniteIndex] <- infiniteACE.info$ACE

    if(method["T1"])
      T1[infiniteIndex] <- infiniteACE.info$T1

    if(method["T2"])
      T2[infiniteIndex] <- infiniteACE.info$T2

    if(hasCustomFun)
      result[infiniteIndex] <- infiniteACE.info$result
  }

  if(withoutCdfs) {
    return(c(if(method["ACE"]) list(ACE=ACE),
             if(method["T1"]) list(T1=T1),
             if(method["T2"]) list(T2=T2),
             if(hasCustomFun) list(result=result)))
  }
  
  if(withoutCi) {
    if(!isSlaveMode && GroupReverse)
      cdfs <- list(Fas0=Fas1, Fas1=Fas0, alphahat=alphahat)
    else
      cdfs <- list(Fas0=Fas0, Fas1=Fas1, alphahat=alphahat)
    
    return(c(if(method["ACE"]) list(ACE=ACE),
             if(method["T1"]) list(T1=T1),
             if(method["T2"]) list(T2=T2),
             if(hasCustomFun) list(result=result),
             cdfs))
  }

  ACE.var.dim <- c(ACE.dim, length(ci.method))
  ACE.var.length <- prod(ACE.var.dim)
  ACE.var.dimnames <- c(list(ACE.dimnames), list(ci.method=ci.method))

  temp <- array(numeric(ACE.var.length), dim=ACE.var.dim,
                                       dimnames=ACE.var.dimnames)
  if(method["ACE"])
    ACE.var <- temp

  if(method['T1'])
    T1.var <- temp

  if(method['T2'])
    T2.var <- temp

  if(hasCustomFun && 'bootstrap' %in% ci.method)
    result.var <- temp[, "bootstrap", drop=FALSE]
  
  ## Done with temp
  rm(temp)
  
  if(!isSlaveMode) {
    ci.map <- vector('list', length(ci.type))
    names(ci.map) <- sprintf("%s", as.character(ci))
  
    for(i in seq_along(ci.type)) {
      if(ci.type[i] == "upper")
        ci.map[[i]] <- ci[i]
      else if(ci.type[i] == "lower")
        ci.map[[i]] <- 1 - ci[i]
      else if(ci.type[i] == "twoSided")
        if(ci[i] < 0.5)
          ci.map[[i]] <- c(ci[i] - ci[i]/2, 1 - ci[i]/2)
        else
          ci.map[[i]] <- c((1-ci[i])/2, ci[i] + (1 - ci[i])/2)
    }

    ci.probs <- sort(unique(unlist(ci.map, recursive=TRUE, use.names=FALSE)))
    ci.probsLen <- length(ci.probs)

    ACE.ci.dim <- c(ACE.dim, ci.probsLen, length(ci.method))
    ACE.ci.length <- prod(ACE.ci.dim)
    ACE.ci.dimnames <- c(list(ACE.dimnames),
                         list(ci.probs= sprintf("%s%%",
                                as.character(ci.probs*100)),
                              ci.method=ci.method))

    ci.map <- lapply(ci.map, ci.probs=ci.probs,
                     ci.prob.names=ACE.ci.dimnames$ci.probs,
                     FUN=function(map, ci.probs, ci.prob.names) {
                       ci.prob.names[match(map, ci.probs)]
                     })
    
    temp <- array(numeric(ACE.ci.length), dim=ACE.ci.dim,
                  dimnames=ACE.ci.dimnames)

    if(method['ACE'])
      ACE.ci <- temp

    if(method['T1'])
      T1.ci <- temp

    if(method['T2'])
      T2.ci <- temp

    if(hasCustomFun && 'bootstrap' %in% ci.method)
      result.ci <- temp[,, "bootstrap", drop=FALSE]
  
    ## Done with temp
    rm(temp)

    ACE.p.dim <- c(ACE.dim, n.test, length(ci.method))
    ACE.p.length <- prod(ACE.p.dim)
    ACE.p.dimnames <- c(list(ACE.dimnames),
                        list(test=names.test, ci.method=ci.method))

    temp <- array(numeric(ACE.p.length), dim=ACE.p.dim,
                  dimnames=ACE.p.dimnames)
    
    if(method['ACE'])
      ACE.p <- temp

    if(method['T1'])
      T1.p <- temp

    if(method['T2'])
      T2.p <- temp

    if(hasCustomFun && 'bootstrap' %in% ci.method)
      result.p <- temp[,, "bootstrap", drop=FALSE]

    ## Done with temp
    rm(temp)
  }
  
  ## run bootstrap method
  if(any(ci.method == 'analytic')) {
    if(method["T1"] || method["T2"])
      stop("analytic ci.method is not compatable with either method T1 or T2")

    if(doInfinite) {
      if(method["ACE"]) {
        ACE.var[infiniteIndex,'analytic'] <- infiniteACE.info$ACE.var[,'analytic']
      }
    }

    if(doFinite) {
      if(method["ACE"]) {
        Omega <- matrix(nrow=5,ncol=5)

        for(i in finiteIndex) {
          .sumCrossUpperTri(Omega) <-
            rbind((p0-s)*(1-z),
                  (p1-s)*z,
                  s*(1/(exp(-beta[i]*y-alphahat[i])+1)-p1/p0)*(1-z),
                  s*(mu0[i]-p0*y/(p1*(exp(-beta[i]*y-alphahat[i])+1)))*(1-z),
                  s*(y1.mean-y)*z)
          
          Omega <- Omega/N
          Omega <- .foldUpperTri(Omega)

          IGamma <-
            solve(matrix(colSums(cbind(1-z,
                                       0,
                                       p1*s*(1-z)/p0^2,
                                       -s*y*(1-z)/(p1*(exp(-beta[i]*y-alphahat[i])+1)),
                                       0,

                                       
                                       0,
                                       z,
                                       -s*(1-z)/p0,
                                       p0*s*y*(1-z)/(p1^2*(exp(-beta[i]*y-alphahat[i])+1)),
                                       0,

                                       
                                       0,
                                       0,
                                       s*exp(-beta[i]*y-alphahat[i])*(1-z)/(exp(-beta[i]*y-alphahat[i])+1)^2,
                                       -p0*s*y*exp(-beta[i]*y-alphahat[i])*(1-z)/(p1*(exp(-beta[i]*y-alphahat[i])+1)^2),
                                       0,

                                       
                                       0,
                                       0,
                                       0,
                                       s*(1-z),
                                       0,

                                       
                                       0,
                                       0,
                                       0,
                                       0,
                                       s*z),
                                 na.rm=TRUE),
                         nrow=5, byrow=FALSE)/N)

          vartheta <- IGamma %*% Omega %*% t(IGamma) / N
          rm(IGamma)
          ACE.var[i,'analytic'] <- vartheta[4,4]+vartheta[5,5] - 2*vartheta[4,5]
        }
        
        ## Done with Omega
        rm(Omega)
      }
    }
    
    if(!isSlaveMode) {
      if(method["ACE"]) {
        calculateCi <- function(i, norm, ACE, sqrt.ACE.var) {
          ACE[i] + norm * sqrt.ACE.var[i]
        }

        ACE.ci[,, 'analytic'] <-
          outer(seq_along(ACE), qnorm(ci.probs),
                FUN=calculateCi, ACE=ACE,
                sqrt.ACE.var=sqrt(ACE.var[, 'analytic']))

        ACE.p[,,'analytic'] <- calc.pvalue(x=ACE, var=ACE.var[,'analytic'],
                                           test=test)
      }
    }
  }

  ## Done with infiniteACE.info

  if(isSlaveMode) {
    cdfs <- list(Fas0=Fas0, Fas1=Fas1, alphahat=alphahat)

    return(c(if(method["ACE"]) list(ACE=ACE, ACE.var=ACE.var),
             if(method["T1"]) list(T1=T1, T1.var=T1.var),
             if(method["T2"]) list(T2=T2, T2.var=T2.var),
             cdfs))
  }

  # Done with T1.var, and T2.var if created
  
  if(any(ci.method == 'bootstrap')) {
    index <- seq_len(N)
    current.fun <- sys.function()
    returned.vals <- c(names.method, if(hasCustomFun) "result")
    Resp.list <- unlist(replicate(N.boot, {
      Resp.vals <- numeric(length(beta)*length(returned.vals))
      dim(Resp.vals) <- c(length(beta), length(returned.vals))
      dimnames(Resp.vals) <- list(beta, returned.vals)
      
      new.index <- sample(N, replace=TRUE)
      new.z <- z[new.index]
      new.s <- s[new.index]
      new.y <- y[new.index]
      if(doInfinite) {
        temp <- sensitivityHHS(z=new.z, s=new.s, y=new.y,
                               bound=bounds,
                               groupings=GroupReverse,
                               ci.method=NULL, custom.FUN=custom.FUN,
                               method=names.method,
                               isSlaveMode=TRUE)[returned.vals]
        for(i in returned.vals) {
          Resp.vals[infiniteIndex,i] <- temp[[i]]
        }
      }

      temp <- current.fun(z=new.z, s=new.s, y=new.y,
                          beta=beta[finiteIndex],
                          groupings=GroupReverse, custom.FUN=custom.FUN,
                          method=names.method,
                          ci.method=NULL, isSlaveMode=TRUE)[returned.vals]

      for(i in returned.vals) {
        Resp.vals[finiteIndex, i] <- temp[[i]]
      }

      Resp.vals
    }, simplify=FALSE), recursive=FALSE)

    ## Done with index and current.fun
    rm(index, current.fun)

    N.bootActual <- length(Resp.list) %/% (length(returned.vals) * ACE.length)

    method.grid <- expand.grid(beta=beta, method=returned.vals)

    ## Done with returned.vals
    rm(returned.vals)
    
    ans <- sapply(split(Resp.list, rep.int(interaction(method.grid),
                                           times=N.bootActual)),
                  function(Resp.vals, probs)
                    c(var(Resp.vals),
                      mean(Resp.vals > 0),
                      mean(Resp.vals < 0),
                      quantile(Resp.vals, probs=probs)),
                    probs=ci.probs)
    ans.ci <- t(ans[c(-1,-2,-3),])
    ans.var <- ans[1,]
    ans.p <- cbind(if(test['upper']) ans[3,],
                   if(test['lower']) ans[2,],
                   if(test['twoSided']) 2 * ifelse(ans[2,] > ans[3,],
                                                   ans[3,], ans[2,]))
    ## Done with ans, Resp.list
    rm(ans, Resp.list)

    if(method["ACE"]) {
      method.index <- method.grid$method == "ACE"
      ACE.ci[,,'bootstrap'] <- ans.ci[method.index,]
      ACE.var[,'bootstrap'] <- ans.var[method.index]
      ACE.p[,,'bootstrap'] <- ans.p[method.index,]
    }

    if(method["T1"]) {
      method.index <- method.grid$method == "T1"
      T1.ci[,,'bootstrap'] <- ans.ci[method.index,]
      T1.p[,,'bootstrap'] <- ans.p[method.index,]
    }

    if(method["T2"]) {
      method.index <- method.grid$method == "T2"
      T2.ci[,,'bootstrap'] <- ans.ci[method.index,]
      T2.p[,,'bootstrap'] <- ans.p[method.index,]
    }

    if(hasCustomFun) {
      method.index <- method.grid$method == "result"
      result.ci[,,'bootstrap'] <- ans.ci[method.index,]
      result.p[,,'bootstrap'] <- ans.p[method.index,]
    }

    ## clean up ans.ci ans.var, and method.grid not used again
    rm(ans.ci, ans.var, ans.p, method.grid, method.index)
  }

  if(!isSlaveMode && GroupReverse)
    cdfs <- list(Fas0=Fas1, Fas1=Fas0[bIndex])
  else
    cdfs <- list(Fas0=Fas0[bIndex], Fas1=Fas1)

  ans <- structure(c(if(method["ACE"]) list(ACE=ACE[bIndex],
                                            ACE.ci=ACE.ci[bIndex,,, drop=FALSE],
                                            ACE.var=ACE.var[bIndex,, drop=FALSE],
                                            ACE.p=ACE.p[bIndex,,, drop=FALSE]),
                     if(method["T1"]) list(T1=T1[bIndex],
                                           T1.ci=T1.ci[bIndex,,, drop=FALSE],
                                           T1.p=T1.p[bIndex,,, drop=FALSE]),
                     if(method["T2"]) list(T2=T2[bIndex],
                                           T2.ci=T2.ci[bIndex,,, drop=FALSE],
                                           T2.p=T2.p[bIndex,,, drop=FALSE]),
                     if(hasCustomFun) list(result=result[bIndex]),
                     if(hasCustomFun && 'bootstrap' %in% ci.method)
                        list(result.ci=result.ci[bIndex,,, drop=FALSE],
                             result.p=result.p[bIndex,,, drop=FALSE]),
                     list(ci.map=ci.map, beta=beta[bIndex],
                          alphahat=alphahat[bIndex]),
                     cdfs),
                   class=c("sensitivity.1.0d", "sensitivity.0d", "sensitivity"),
                   parameters=list(z0=groupings[1], z1=groupings[2],
                     selected=selection, s0=empty.principal.stratum[1],
                     s1=empty.principal.stratum[2]))
  if('bootstrap' %in% ci.method)
    attr(ans, 'N.boot') <- N.boot
  
  return(ans)
}
