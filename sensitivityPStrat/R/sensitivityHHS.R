.PrepDataObj <- function(z, s, y, GroupReverse, withoutCdfs, method,
                         custom.FUN) {
  ## Various summary numbers
  ## N : total number of records in data
  ## N0, N1           : number of records in the first and second groups
  ## n0, n1           : number of selected records in the first and second group
  ## p0, p1           : probablity of a record being selected in the frist and
  ##                       second group
  ## z0.s1, z1.s0     : index vector of the selected records in the first and
  ##                       second groups
  ## y0, y1           : y values of selected records of the first and second
  ##                       groups
  ## RR               : Ratio of probablities
  ## VE               : Vacine Efficecy
  ## R                : Rank of y of selected records
  ## Ras0             : distribution of the Rank of the first group
  ## F0, F1           : distribution of the first and second groups
  ## y0.uniq, y1.uniq : the sorted unique values of y in the selected records
  ##                      of the first and second groups

  z0.s1 <- !z & s
  z1.s1 <- z & s

  y0 <- y[z0.s1]
  y1 <- y[z1.s1]

  N <- length(z)
  N0 <- sum(!z)
  n0 <- sum(z0.s1)
  p0 <- n0/N0

  N1 <- sum(z)
  n1 <- sum(z1.s1)
  p1 <- n1/N1

  RR <- min(p1/p0, 1L)

  VE <- 1L - RR

  R <- rank(y[s])

  Fn0 <- ecdf(y0)
  Fn1 <- ecdf(y1)

  y0.mean <- mean(y0)
  y1.mean <- mean(y1)

  c(list(withoutCdfs=withoutCdfs, GroupReverse=GroupReverse,
         z=z, s=s, y=y, z0.s1=z0.s1, z1.s1=z1.s1, N=N, R1=R[z[s]],
         N0=N0, N1=N1, n0=n0, n1=n1, p0=p0, p1=p1, RR=RR, VE=VE, Fn0=Fn0,
         Fn1=Fn1, y0=y0, y1=y1, method=method, R0.uniq=sort(unique(R[!z[s]])),
         y0.mean=y0.mean, y1.mean=y1.mean, custom.FUN=custom.FUN),
    eval(expression(list(y0.uniq=x, F0=y)), envir=environment(Fn0)),
    eval(expression(list(y1.uniq=x, F1=y)), envir=environment(Fn1)))
}

.CalcDeltaMu <- function(mu0, mu1, GroupReverse) {
  if(GroupReverse) {
    mu0 - mu1
  } else {
    mu1 - mu0
  }
}

.CreateACERetValHHS <- function(Fas0, obj) {
  dFas0 <- diff(c(0L, Fas0))

  if(obj$method['ACE'] || obj$method['T2']) {
    mu0 <- sum(obj$y0.uniq * dFas0)
  }
  
  if(obj$method['ACE']) {
    ACE <- .CalcDeltaMu(mu0=mu0, mu1=obj$y1.mean, GroupReverse=obj$GroupReverse)
  }

  if(obj$method['T1']) {
    Rmu0 <- sum(obj$R0.uniq * dFas0)
    Rmu1 <- mean(obj$R1)

    T1 <- .CalcDeltaMu(mu0=Rmu0, mu1=Rmu1, GroupReverse=obj$GroupReverse)
  }

  if(obj$method['T2']) {
    indx <- rep(c(FALSE, TRUE), times=c(obj$n0, obj$n1))
    ystar0 <- obj$y0 - obj$y0.mean - mu0
    ystar1 <- obj$y1

    Rstar <- rank(c(ystar0,ystar1))

    T2 <- .CalcDeltaMu(mu0=mean(Rstar[!indx]), mu1=mean(Rstar[indx]),
                       GroupReverse=obj$GroupReverse)
  }

  FnAs0 <- stepfun(x=obj$y0.uniq, y=c(0, Fas0))
  
  if(!is.null(obj$custom.FUN)) {
    result <- obj$custom.FUN(mu0=mu0, mu1=obj$y1.mean,
                             p0=obj$p0, p1=obj$p1)
  }
  
  resp <- c(if(obj$method['ACE']) list(ACE=ACE),
            if(obj$method['T1']) list(T1=T1),
            if(obj$method['T2']) list(T2=T2),
            if(!is.null(obj$custom.FUN)) list(result=result))
  
  if(obj$withoutCdfs)
    return(resp)
  
  return(c(resp,
           list(Fas0=Fas0, FnAs0=FnAs0)))
}

.CalcUpperACEHHS <- function(obj) {
  ## Variables
  ## qc : y0 value that represents the RRth percentila of y0
  ## Fas0, Fas1 : distribution funtion for the always selected stratum for the first and second groups
  ## dFas0, dFas1 : delta of Fas0 and Fas1

  qc <- quantile(obj$y0, probs=obj$RR)
  Fas0 <- ifelse(obj$y0.uniq <= qc & obj$F0 < obj$RR, obj$F0/obj$RR, 1)
  
  .CreateACERetValHHS(Fas0=Fas0, obj=obj)
}
  
.CalcLowerACEHHS <- function(obj) {
  ## Variables
  ## qc : y0 value that represents the RRth percentila of y0
  ## Fas0, Fas1 : distribution funtion for the always selected stratum for the first and second groups
  ## dFas0, dFas1 : delta of Fas0 and Fas1

  qc <- quantile(obj$y0, probs=obj$VE)
  Fas0 <- ifelse(obj$y0.uniq >= qc, (obj$F0 - obj$VE)/obj$RR, 0)
  
  .CreateACERetValHHS(Fas0=Fas0, obj=obj)
}  

sensitivityHHS <- function(z, s, y, bound=c("upper","lower"),
                           selection, groupings, empty.principal.stratum,
                           ci=0.95, ci.method=c("bootstrap", "analytic"),
                           ci.type="twoSided", custom.FUN=NULL,
                           na.rm=FALSE, N.boot=100, upperTest=FALSE,
                           lowerTest=FALSE, twoSidedTest=TRUE,
                           method=c("ACE", "T1", "T2"),
                           isSlaveMode=FALSE)
{
  withoutCdfs <- isSlaveMode && !missing(ci.method) && is.null(ci.method)
  withoutCi <- ((!isSlaveMode && !missing(ci.method) && ci.method == "") ||
                (isSlaveMode && !(!missing(ci.method) &&
                                 !is.null(ci.method) &&
                                 'analytic' %in% ci.method)))
  
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

    if(any(is.na(z) | is.na(s)))
      stop("s, z cannot contain any NA values")
    
    if(any(s & is.na(y)))
      stop("selected y values cannot be NA")

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
  
  bound.orig <- bound
  bound <- unique(match.arg(bound, several.ok=TRUE))
  boundIndex <- match(bound.orig, bound)

  hasCustomFun <- !is.null(custom.FUN)
  
  if(withoutCi)
    ci.method <- NULL
  else if(isSlaveMode)
    ci.method <- "analytic"
  else if(missing(ci.method) || is.null(ci.method)) {
    ci.method <- 'bootstrap'
  } else {
    ci.method <- sort(unique(match.arg(ci.method, several.ok=TRUE)))
  }

  method <- match.bitarg(method)
  n.method <- sum(method)
  names.method <- names(method)[method]

  test <- c(upper=upperTest, lower=lowerTest, twoSided=twoSidedTest)
  n.test <- sum(test)
  names.test <- names(test)[test]
  
  UpperIndex <- bound == "upper"
  LowerIndex <- bound == "lower"

  DoUpper <- any(UpperIndex)
  DoLower <- any(LowerIndex)
  
  datObj <- .PrepDataObj(z=z, s=s, y=y, GroupReverse=GroupReverse,
                         withoutCdfs=withoutCdfs, method=method,
                         custom.FUN=custom.FUN)

  ACE.dim <- length(bound)
  ACE.length <- prod(ACE.dim)
  ACE.dimnames <- bound
  
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

  ## Done with temp
  rm(temp)
  
  if(!withoutCdfs) {
    FnAs0 <- funVector(ACE.length)
    names(FnAs0) <- ACE.dimnames

    FnAs1 <- funVector(1L)
    FnAs1[1] <- datObj$Fn1
  }
  
  if(DoUpper) {
    UpperObj <- .CalcUpperACEHHS(datObj)

    if(method["ACE"])
      ACE['upper'] <- UpperObj$ACE

    if(method['T1'])
      T1['upper'] <- UpperObj$T1

    if(method['T2'])
      T2['upper'] <- UpperObj$T2

    if(hasCustomFun)
      result['upper'] <- UpperObj$result

    if(!withoutCdfs) {
      FnAs0['upper'] <- UpperObj$FnAs0
    }

    ## Done with UpperObj
    rm(UpperObj)
  }

  if(DoLower) {
    LowerObj <- .CalcLowerACEHHS(datObj)

    if(method["ACE"])
      ACE['lower'] <- LowerObj$ACE
    
    if(method['T1'])
      T1['lower'] <- LowerObj$T1

    if(method['T2'])
      T2['lower'] <- LowerObj$T2

    if(hasCustomFun)
      result['lower'] <- LowerObj$result

    if(!withoutCdfs)
      FnAs0['lower'] <- LowerObj$FnAs0

    ## Done with LowerObj
    rm(LowerObj)
  }

  if(withoutCdfs)
    return(c(if(method["ACE"]) list(ACE=ACE),
             if(method["T1"]) list(T1=T1),
             if(method["T2"]) list(T2=T2),
             if(hasCustomFun) list(result=result)))

  if(!isSlaveMode && GroupReverse) {
    Fas0 <- FnAs1
    Fas1 <- FnAs0
  } else {
    Fas0 <- FnAs0
    Fas1 <- FnAs1
  }
  
  if(withoutCi) {
    return(c(if(method["ACE"]) list(ACE=ACE),
             if(method["T1"]) list(T1=T1),
             if(method["T2"]) list(T2=T2),
             if(hasCustomFun) list(result=result),
             list(Fas0=Fas0, Fas1=Fas1)))
  }

  ACE.var.dim <- c(ACE.dim, length(ci.method))
  ACE.var.length <- prod(ACE.var.dim)
  ACE.var.dimnames <- c(list(ACE.dimnames), list(ci.method=ci.method))

  temp <- array(numeric(ACE.var.length), dim=ACE.var.dim,
                dimnames=ACE.var.dimnames)

  if(method["ACE"])
    ACE.var <- temp
  
  if(method["T1"])
    T1.var <- temp

  if(method["T2"])
    T2.var <- temp

  if(hasCustomFun && 'bootstrap' %in% ci.method) {
    result.var <- temp[, 'bootstrap', drop=FALSE]
  }
  
  ## Done with temp
  rm(temp)

  ## Do analytic method
  if('analytic' %in% ci.method) {
    .FeatureNotYetImplemented("analytic method")
  }

  if(isSlaveMode) {
    return(c(if(method["ACE"]) list(ACE=ACE, ACE.var=ACE.var),
             if(method["T1"]) list(T1=T1, T1.var=T1.var),
             if(method["T2"]) list(T2=T2, T2.var=T2.var),
             if(method["result"]) list(result=result),
             list(Fas0=Fas0, Fas1=Fas1)))
  }

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

  if(hasCustomFun && 'bootstrap' %in% ci.method) {
    result.p <- temp[,, 'bootstrap', drop=FALSE]
  }

  # Done with temp
  rm(temp)
  
  ci.map <- vector('list', length(ci.type))
  names(ci.map) <- ci
  
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
  ACE.ci.dimnames <- c(list(bound=ACE.dimnames),
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

  if(method["ACE"])
    ACE.ci <- temp

  if(method['T1'])
    T1.ci <- temp
  
  if(method['T2'])
    T2.ci <- temp

  if(hasCustomFun && ci.method %in% 'bootstrap')
    result.ci <- temp[,,'bootstrap', drop=FALSE]
  
  ## Done with temp
  rm(temp)
  
  if('analytic' %in% ci.method) {
    calculateCi <- function(i, norm, ACE, sqrt.ACE.var) {
      ACE[i] + norm * sqrt.ACE.var[i]
    }


    if(method["ACE"]) {
      ACE.ci[,,'analytic'] <- outer(seq_along(ACE), qnorm(ci.probs),
                                    FUN=calculateCi, ACE=ACE,
                                    sqrt.ACE.var=sqrt(ACE.var[,'analytic']))
      ACE.p[,,'analytic'] <- calc.pvalue(x=ACE, var=ACE.var[,'analytic'],
                                         test=test)
    }
  }

  ## run bootstrap method
  if(any(ci.method == 'bootstrap')) {
    bootACECalc <- function(i, nVal, z, s, y, bound, GroupReverse,
                            names.method, returned.vals, current.fun,
                            custom.FUN) {
      index <- sample(nVal, replace=TRUE)

      ans <- current.fun(z=z[index], s=s[index], y=y[index], bound=bound,
                         groupings=GroupReverse, ci.method=NULL,
                         custom.FUN=custom.FUN,
                         method=names.method, isSlaveMode=TRUE)[returned.vals]

      return(ans)
    }

    returned.vals <- c(names.method, if(hasCustomFun) "result")
    Resp.list <- unlist(lapply(integer(N.boot), nVal=datObj$N,
                               z=datObj$z, s=datObj$s, y=datObj$y, bound=bound,
                               GroupReverse=GroupReverse,
                               names.method=names.method,
                               returned.vals=returned.vals,
                               custom.FUN=custom.FUN,
                               current.fun=sys.function(), FUN=bootACECalc))

    N.bootActual <- length(Resp.list) %/% (length(returned.vals) * ACE.length)
    method.grid <- expand.grid(bound=bound, method=returned.vals)
    ans <- sapply(split(Resp.list, rep.int(interaction(method.grid),
                                           times=N.bootActual)),
                  FUN=function(Resp.vals, probs)
                    c(var(Resp.vals),
                      mean(Resp.vals > 0),
                      mean(Resp.vals < 0),
                      quantile(Resp.vals, probs=probs)),
                  probs=ci.probs)
    ans.ci <- t(ans[c(-1,-2,-3),, drop=FALSE])
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

    ## clean up ans.ci ans.var, ans.p, and method.grid not used again
    rm(ans.ci, ans.var, ans.p, method.grid)
  }

  ans <-
    structure(c(if(method["ACE"]) list(ACE=ACE[boundIndex],
                                       ACE.ci=ACE.ci[boundIndex,,, drop=FALSE],
                                       ACE.var=ACE.var[boundIndex,, drop=FALSE],
                                       ACE.p=ACE.p[boundIndex,,, drop=FALSE]),
                if(method["T1"]) list(T1=T1[boundIndex],
                                      T1.ci=T1.ci[boundIndex,,, drop=FALSE],
                                      T1.p=T1.p[boundIndex,,, drop=FALSE]),
                if(method["T2"]) list(T2=T2[boundIndex],
                                      T2.ci=T2.ci[boundIndex,,, drop=FALSE],
                                      T2.p=T2.p[boundIndex,,, drop=FALSE]),
                if(hasCustomFun) list(result=result[boundIndex],
                                      result.ci=result.ci[boundIndex,,, drop=FALSE],
                                      result.p=result.p[boundIndex,,, drop=FALSE]),
                list(ci.map=ci.map, bound=bound[boundIndex],
                     beta=c(-Inf, Inf)[match(bound, c('upper','lower'))][boundIndex],
                     Fas0=Fas0[boundIndex], Fas1=Fas1[boundIndex])),
              class=c("sensitivity.1.0d", "sensitivity.0d", "sensitivity"),
              parameters=list(z0=groupings[1], z1=groupings[2],
                selected=selection, s0=empty.principal.stratum[1],
                s1=empty.principal.stratum[2]))

  if('bootstrap' %in% ci.method) {
    attr(ans, 'N.boot') <- N.boot
    attr(ans, 'N.bootActual') <- N.bootActual
  }

  return(ans)
}
