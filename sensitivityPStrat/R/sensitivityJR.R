sensitivityJR <- function(z, s, y, beta0, beta1, phi, Pi, psi,
                          selection, groupings, ci=0.95,
                          ci.method=c("analytic", "bootstrap"),
                          ci.type="twoSided", custom.FUN=NULL, na.rm=FALSE,
                          N.boot=100, interval=c(-100,100),
                          upperTest=FALSE, lowerTest=FALSE, twoSidedTest=TRUE,
                          verbose=getOption("verbose"), isSlaveMode=FALSE)
{
  withoutCdfs <- isSlaveMode && !missing(ci.method) && is.null(ci.method)
  withoutCi <- ((!isSlaveMode && !missing(ci.method) && ci.method == "") ||
                (isSlaveMode && !(!missing(ci.method) &&
                                 !is.null(ci.method) &&
                                 'analytic' %in% ci.method)))
  hasCustomFun <- !is.null(custom.FUN)

  calc.coefs <- function(y, beta, dF, RR, interval) {
    coefs <- vector(length(beta), "list")

    for(i in seq_along(beta)) {
      alphahat <- .calc.alphahat(beta.y=beta[i]*y, dF=dF, C=RR,
                                 interval=interval)
      w <- .calc.w(alpha=alphahat, beta.y=beta[i]*y)
      coefs[[i]] <- list(alphahat=alphahat, w=w)
    }

    return(coefs)
  }

  if(!isSlaveMode) {
    ## Not running a boot strap mode
    ## Run error checks on variables.
    ErrMsg <- NULL
    ErrMsg <- c(.CheckSelection(selection, s),
                .CheckGroupings(groupings),
                .CheckPhiPiPsi(phi=phi, Pi=Pi, psi=psi),
                .CheckLength(z=z, s=s, y=y),
                .CheckZ(z, groupings, na.rm=na.rm),
                .CheckS(s, na.rm=na.rm),
                .CheckY(y, s, selection, na.rm=na.rm),
                .CheckCi(ci=ci, ci.type=ci.type))
    
    if(length(ErrMsg) > 0L)
      stop(paste(ErrMsg, collapse="\n  "))

    s <- s == selection
    
    z <- z == groupings[2L]

    if(na.rm == TRUE) {
      naIndex <- !(is.na(s) | is.na(z) | (s & (is.na(y))))

      z <- z[naIndex]
      s <- s[naIndex]
      y <- y[naIndex]
    }

    if(any(is.na(z) | is.na(s)))
      stop("s, z cannot contain any NA values")
    
    if(any(s & is.na(y)))
      stop("selected y values cannot be NA")

    ErrMsg <- .CheckPhiPiPsi(phi=phi, Pi=Pi, psi=psi, p0=sum(!z & s)/sum(!z), p1=sum(z & s)/sum(z))

    if(length(ErrMsg) > 0L)
      stop(paste(ErrMsg, collapse="\n  "))

    if(missing(ci.type)) {
      ci.type <- rep('twoSided', length.out=length(ci))
    } else {
      ci.type <- match.arg(ci.type, c('upper', 'lower', 'twoSided'),
                           several.ok=TRUE)
    }
  }

  test <- c(upper=upperTest, lower=lowerTest, twoSided=twoSidedTest)
  n.test <- sum(test)
  names.test <- names(test)[test]
  
  if(withoutCi)
    ci.method <- NULL
  else if(isSlaveMode)
    ci.method <- "analytic"
  else
    ci.method <- sort(unique(match.arg(ci.method, several.ok=TRUE)))

  z0.s1 <- !z & s
  z1.s1 <- z & s
  
  y0 <- y[z0.s1]
  y1 <- y[z1.s1]

  y0.mean <- mean(y0)
  y1.mean <- mean(y1)

  N <- length(z)
  N0 <- sum(!z)
  n0 <- sum(z0.s1)
  p0 <- n0/N0

  N1 <- sum(z)
  n1 <- sum(z1.s1)
  p1 <- n1/N1

  RR <- p1/p0
  VE <- 1 - RR

  ## Calc Pi because it will be needed later
  tmp <- .calcPiPhiPsi(Pi=Pi, phi=phi, psi=psi, p0=p0, p1=p1)
  Pi <- tmp$Pi
  psi <- tmp$psi
  phi <- tmp$phi
  sens.var <- tmp$sens.var
  

  Fn0 <- ecdf(y0)
  y0.uniq <- knots(Fn0)
  F0 <- Fn0(y0.uniq)
  dF0 <- diff(c(0, F0))
  
  Fn1 <- ecdf(y1)
  y1.uniq <- knots(Fn1)
  F1 <- Fn1(y1.uniq)
  dF1 <- diff(c(0, F1))

  ACE.dim <- c(length(beta0), length(beta1), length(Pi))
  ACE.length <- prod(ACE.dim)
  ACE.dimnames <- list(format(beta0, trim=TRUE),
                       format(beta1, trim=TRUE),
                       format(switch(sens.var,
                                     Pi=Pi,
                                     phi=phi,
                                     psi=psi),
                              trim=TRUE, digits=4,
                              drop0trailing=TRUE))
  names(ACE.dimnames) <- c("beta0", "beta1", sens.var)

  ACE <- array(numeric(ACE.length), dim=ACE.dim, dimnames=ACE.dimnames)

  if(hasCustomFun)
    result <- ACE
  
  FnAs0.dim <- ACE.dim[-2L]
  FnAs0.length <- prod(FnAs0.dim)
  FnAs0.dimnames <- ACE.dimnames[-2L]

  mu0 <- alphahat0 <- array(numeric(FnAs0.length), dim=FnAs0.dim,
                            dimnames=FnAs0.dimnames)

  FnAs1.dim <- ACE.dim[-1L]
  FnAs1.length <- prod(FnAs1.dim)
  FnAs1.dimnames <- ACE.dimnames[-1L]

  mu1 <- alphahat1 <- array(numeric(FnAs1.length), dim=FnAs1.dim,
                            dimnames=FnAs1.dimnames)

  if(!isSlaveMode) {
    FnAs0 <- funArray(dim=FnAs0.dim,
                      dimnames=FnAs0.dimnames)

    FnAs1 <- funArray(dim=FnAs1.dim,
                      dimnames=FnAs1.dimnames)
  }

  a0 <- Pi/p0
  a1 <- Pi/p1
    
  for(i in seq_along(Pi)) {
    if(phi[i] == 1) {
      ACE.info <- sensitivityGBH(z=z,s=s,y=y,beta=beta0,
                                 groupings=FALSE,
                                 ci.method=ci.method,
                                 custom.FUN=custom.FUN, isSlaveMode=TRUE)
      
      if(!withoutCdfs) {
        alphahat0[,i] <- ACE.info$alphahat
        FnAs0[,i] <- ACE.info$FnAs0
        alphahat1[,i] <- NA
        FnAs1[,i] <- ACE.info$FnAs1
      }        
      
      ACE[,,i] <- ACE.info$ACE

      if(hasCustomFun)
        result[,,i] <- ACE.info$result
      
      next
    }
    
    if(Pi[i] == 0) {
      ACE[,,i] <- NA
      alphahat0[,i] <- NA
      alphahat1[,i] <- NA

      if(hasCustomFun)
        result[,,i] <- NA
      
      next
    }
    
    q0c <- quantile(y0, probs=c(a0[i], 1-a0[i]))
    q1c <- quantile(y1, probs=c(a1[i], 1-a1[i]))
  
    for(j in seq_along(beta0)) {
      if(is.finite(beta0[j])) {
        alphahat0[j,i] <- .calc.alphahat(beta.y=beta0[j]*y0.uniq, dF=dF0,
                                         C=a0[i],
                                         interval=interval)
        w0 <- .calc.w(alpha=alphahat0[j,i], beta.y=beta0[j]*y0.uniq)

        dFas0 <- dF0*w0/a0[i]
        if(!isSlaveMode)
          Fas0 <- cumsum(dFas0)
        
        mu0[j,i] <- sum(y0.uniq * dFas0)
      } else {
        if(beta0[j] == Inf) {
          Fas0 <- ifelse(y0.uniq <= q0c[1L] & F0 < a0[i], F0/a0[1], 1)
        } else if(beta0[j] == -Inf) {
          Fas0 <- ifelse(y0.uniq >= q0c[2L], (F0 - (1 - a0[i]))/a0[i], 0)
        } else {
          stop("Invalid beta0 value ", beta0[j])
        }

        alphahat0[j,i] <- NA
        mu0[j,i] <- sum(y0.uniq * diff(c(0, Fas0)))
      }

      if(!isSlaveMode) {
        FnAs0[j,i] <- stepfun(y0.uniq, c(0, Fas0))
      }
    }
    
    for(j in seq_along(beta1)) {
      if(is.finite(beta1[j])) {
        alphahat1[j,i] <- .calc.alphahat(beta.y=beta1[j]*y1.uniq, dF=dF1,
                                         C=a1[i],
                                         interval=interval)
        w1 <- .calc.w(alpha=alphahat1[j,i], beta.y=beta1[j]*y1.uniq)

        dFas1 <- dF1*w1/a1[i]

        if(!isSlaveMode)
          Fas1 <- cumsum(dFas1)
        
        mu1[j,i] <- sum(y1.uniq * dFas1)
      } else if(is.infinite(beta1[j])) {
        if(beta1[j] == Inf) {
          Fas1 <- ifelse(y1.uniq <= q1c[1L] & F1 < a1[i], F1/a1[i], 1L)
        } else if(beta1[j] == -Inf) {
          Fas1 <- ifelse(y1.uniq >= q1c[2L], (F1 - a1[i])/(a1[i]), 0L)
        }

        alphahat1[j,i] <- NA
        mu1[j,i] <- sum(y1.uniq * diff(c(0L, Fas1)))
      } else {
        stop("Invalid beta1 value ", beta1[j])
      }

      if(!isSlaveMode)
        FnAs1[j,i] <- stepfun(y1.uniq, c(0, Fas1))
    }

    if(hasCustomFun) {
      for(j in seq_along(beta0)) {
        for(k in seq_along(beta1)) {
          result[j,k,i] <- custom.FUN(mu0=mu0[j,i], mu1=mu1[k,i], p0=p0, p1=p1)
        }
      }
    }
    
    ACE[,,i] <- outer(mu0[,i], mu1[,i], function(mu0, mu1) mu1 - mu0)
  }
  
  if(withoutCdfs)
    return(c(list(ACE=ACE),
             if(hasCustomFun) list(result=result)))

  cdfs<-list(beta0=beta0, alphahat0=alphahat0, Fas0=FnAs0,
             beta1=beta1, alphahat1=alphahat1, Fas1=FnAs1,
             phi=phi, Pi=Pi, psi=psi)

  if(withoutCi)
    return(c(list(ACE=ACE),
             if(hasCustomFun) list(result=result),
             cdfs))
  
  if(!isSlaveMode) {
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
    ACE.ci.dimnames <- c(ACE.dimnames,
                         list(ci.prob=sprintf("%s%%",
                                as.character(ci.probs*100)),
                              ci.method=ci.method))
  
    ci.map <- lapply(ci.map, ci.probs=ci.probs,
                     ci.prob.names=ACE.ci.dimnames$ci.probs,
                     FUN=function(map, ci.probs, ci.prob.names) {
                       ci.prob.names[match(map, ci.probs)]
                     })
    
    ACE.ci <- array(numeric(ACE.ci.length), dim=ACE.ci.dim,
                    dimnames=ACE.ci.dimnames)

    if(hasCustomFun && 'bootstrap' %in% ci.method) {
      result.ci <- ACE.ci[,,,,"bootstrap", drop=FALSE]
    }
  }
  
  ACE.var.dim <- c(ACE.dim, length(ci.method))
  ACE.var.length <- prod(ACE.var.dim)
  ACE.var.dimnames <- c(ACE.dimnames, list(ci.method=ci.method))

  ACE.var <- array(numeric(ACE.var.length), dim=ACE.var.dim,
                   dimnames=ACE.var.dimnames)

  ACE.p.dim <- c(ACE.dim, n.test, length(ci.method))
  ACE.p.length <- prod(ACE.p.dim)
  ACE.p.dimnames <- c(ACE.dimnames, list(test=names.test, ci.method=ci.method))
  
  ACE.p <- array(numeric(ACE.p.length), dim=ACE.p.dim, dimnames=ACE.p.dimnames)

  if(hasCustomFun && 'bootstrap' %in% ci.method) {
    result.var <- ACE.var[,,,"bootstrap", drop=FALSE]
    result.p <- ACE.p[,,,,"bootstrap", drop=FALSE]
  }
  
  ## run bootstrap method
  if(any(ci.method == "analytic")) {
    if(verbose)
      cat("running Analytic")
    
    Omega <- matrix(nrow=6,ncol=6)
    for(k in seq_along(Pi)) {
      if(phi[k] == 1) {
        ACE.var[,, k, "analytic"] <- ACE.info$ACE.var
      } else for(i in seq_along(beta0)) {
        for(j in seq_along(beta1)) {
          U <- rbind((1-z)*(p0 - s),
                     (z)*(p1 - s),
                     (1-z)*s*(1/(1+exp(-alphahat0[i,k] - beta0[i]*y)) - Pi[k]/p0),
                     (z)*s*(1/(1+exp(-alphahat1[j,k] - beta1[j]*y)) - Pi[k]/p1),
                     (1-z)*s*(mu0[i,k] - y*p0/Pi[k]/(1+exp(-alphahat0[i,k] - beta0[i]*y))),
                     (z)*s*(mu1[j,k] - y*p1/Pi[k]/(1+exp(-alphahat1[j,k] - beta1[j]*y))))
          .sumCrossUpperTri(Omega) <- U
          Omega <- Omega / N
          Omega[lower.tri(Omega)] <- t(Omega)[lower.tri(Omega)]

          Gamma <- matrix(colSums(cbind(1-z,
                                        0,
                                        s*(1-z)*Pi[k]/p0^2,
                                        0,
                                        -s*y*(1-z)/Pi[k]/(exp(-beta0[i]*y-alphahat0[i,k])+1),
                                        0,

                                        
                                        0,
                                        z,
                                        0,
                                        s*z*Pi[k]/p1^2,
                                        0,
                                        -s*y*z/Pi[k]/(exp(-beta1[j]*y-alphahat1[j,k])+1),

                                        
                                        0,
                                        0,
                                        s*exp(-beta0[i]*y-alphahat0[i,k])*(1-z)/(exp(-beta0[i]*y-alphahat0[i,k])+1)^2,
                                        0,
                                        -p0*s*y*exp(-beta0[i]*y-alphahat0[i,k])*(1-z)/Pi[k]/(exp(-beta0[i]*y-alphahat0[i,k])+1)^2,
                                        0,

                                        
                                        0,
                                        0,
                                        0,
                                        s*exp(-beta1[j]*y-alphahat1[j,k])*z/(exp(-beta1[j]*y-alphahat1[j,k])+1)^2,
                                        0,
                                        -p1*s*y*exp(-beta1[j]*y-alphahat1[j,k])*z/Pi[k]/(exp(-beta1[j]*y-alphahat1[j,k])+1)^2,

                                        
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
                                        0,
                                        s*z), na.rm=TRUE),
                          nrow=6,byrow=TRUE) / N

          IGamma <- solve(Gamma)
          ## vartheta <- tcrossprod(IGamma %*% Omega, IGamma) / N
          vartheta <- (t(IGamma) %*% Omega %*% IGamma) / N
          ## cat(all.equal(vartheta, alt), '\n',sep='')
          ACE.var[i, j, k, "analytic"] <- vartheta[5, 5] + vartheta[6, 6] - 2 * vartheta[5, 6]

          if(verbose) cat(".")
        }
      }
    }

    calculateCi <- function(i, norm, ACE, sqrt.ACE.var) {
      ACE[i] + norm * sqrt.ACE.var[i]
    }

    if(!isSlaveMode) {
      ACE.ci[,,,,"analytic"] <- outer(seq_along(ACE), qnorm(ci.probs),
                                      FUN=calculateCi, ACE=ACE,
                                      sqrt.ACE.var=sqrt(ACE.var[,,,'analytic']))
      ACE.p[,,,, 'analytic'] <- calc.pvalue(x=ACE, var=ACE.var[,,,'analytic'],
                                            test=test)
    }
    
    if(verbose) cat("\n")
  }

  if(isSlaveMode) {
    return(c(list(ACE=ACE, ACE.var=ACE.var), cdfs))
  }
    
  if(any(ci.method == "bootstrap")) {
    if(verbose) cat("running Bootstrap")
    ACE.list <- array(dim=c(ACE.dim,N.boot))
    if(hasCustomFun)
      result.list <- ACE.list
    
    for(i in seq_len(N.boot)) {
      new.index <- sample(seq_len(N), N, replace=TRUE)
      ans <- Recall(z=z[new.index],s=s[new.index],y=y[new.index],
                    beta0=beta0, beta1=beta1, psi=psi,
                    groupings=FALSE, ci.method=NULL,
                    custom.FUN=custom.FUN, isSlaveMode=TRUE)
      ACE.list[,,,i] <- ans$ACE

      if(hasCustomFun)
        result.list[,,,i] <- ans$result
      
      if(verbose) cat(".")
    }

    for(k in seq_along(Pi)) {
      for(i in seq_along(beta0)) {
        for(j in seq_along(beta1)) {
          ACE.ci[i,j,k,,'bootstrap'] <- quantile(ACE.list[i,j,k,], probs=ci.probs)
          ACE.var[i,j,k,'bootstrap'] <- var(ACE.list[i,j,k,])
          lower <- mean(ACE.list[i,j,k,] > 0)
          upper <- mean(ACE.list[i,j,k,] < 0)

          if(upperTest)
            ACE.p[i,j,k,'upper','bootstrap'] <- upper

          if(lowerTest)
            ACE.p[i,j,k,'lower','bootstrap'] <- lower

          if(twoSidedTest)
            ACE.p[i,j,k,'twoSided', 'bootstrap'] <- 2 * min(lower, upper)

          if(hasCustomFun) {
            result.ci[i,j,k,,'bootstrap'] <- quantile(result.list[i,j,k,], probs=ci.probs)
            result.var[i,j,k,'bootstrap'] <- var(result.list[i,j,k,])
            lower <- mean(result.list[i,j,k,] > 0)
            upper <- mean(result.list[i,j,k,] < 0)

            if(upperTest)
              result.p[i,j,k,'upper','bootstrap'] <- upper

            if(lowerTest)
              result.p[i,j,k,'lower','bootstrap'] <- lower

            if(twoSidedTest)
              result.p[i,j,k,'twoSided', 'bootstrap'] <- 2 * min(lower, upper)
          }
        }
      }
    }
    if(verbose) cat("\n")
  }

  return(structure(c(list(ACE=ACE, ACE.ci=ACE.ci, ACE.var=ACE.var, ACE.p=ACE.p),
                     if(hasCustomFun) list(result=result),
                     if(hasCustomFun && 'bootstrap' %in% ci.method)
                         list(result.ci=result.ci, result.var=result.var,
                              result.p=result.p),
                     list(ci.map=ci.map),
                     cdfs),
                   class=c("sensitivity.2.0d", "sensitivity.0d", "sensitivity"),
                   N.boot=N.boot,
                   parameters=list(z0=groupings[1], z1=groupings[2],
                     selected=selection)))
}
