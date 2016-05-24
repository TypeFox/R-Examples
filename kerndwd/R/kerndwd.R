kerndwd = function(x, y, kern, family=c("DWD", "LUM", "other"),
  lambda, qval=1, aval=1, cval=1, const, dev, wt, 
  eps=1e-05, maxit=1e+05) {
  ####################################################################
  ## data setup
  this.call = match.call()
  if (length(levels(factor(y))) == 2)
    y = c(-1, 1)[as.factor(drop(y))]
  x = as.matrix(x)
  if (!all(y %in% c(-1, 1))) 
    stop("y should be a factor with two levels.")
  nobs = as.integer(NROW(x))
  np = as.integer(NCOL(x))
  if (length(y) != nobs) 
    stop("x and y have different number of observations.")
  if (missing(kern)) {
    kern = rbfdot(sigma=1)
    cat("'kern' is missing: Gaussian kernel is used.\n")
  }
  ####################################################################
  # qval = as.double(qval)
  maxit = as.integer(maxit)
  eps = as.double(eps)
  gam = as.double(1e-8) # a tiny value to avoid matrix sigularity
  ####################################################################
  ## lambda setup
  if (missing(lambda)) {
    stop("Users have to provide a lambda sequence.")
  } else {
    ulam = as.double(rev(sort(lambda)))
    nlam = as.integer(length(lambda))
  }
  family = match.arg(family)
  if (!match(family, c("DWD", "LUM", "other"), FALSE)) {
    warning("Only 'DWD', 'LUM', and 'other' are available; 
      'DWD' used.")
    family = "DWD"
  } 
  fit = switch(family,
      DWD = dwdpath(x, y, nobs, np, kern, 
        qval, ulam, nlam, wt, eps, maxit, gam),
      LUM = lumpath(x, y, nobs, np, kern, 
        aval, cval, ulam, nlam, wt, eps, maxit, gam),
      other = otherpath(x, y, nobs, np, kern, 
      const, dev, ulam, nlam, wt, eps, maxit, gam) 
  )
  fit.call = this.call
  class(fit) = c(class(fit), "kerndwd")
  fit
} 

dwdpath = function(x, y, nobs, np, kern, qval, 
  ulam, nlam, wt, eps, maxit, gam) {
  #################################################################### 
  ## check index 
  if (qval <= 0) {
    warning("qval must be positive; set to 1.")
    qval = 1 
  }  
  ####################################################################
  if (missing(wt) || is.null(wt)) {
    if (class(kern)[[1]] == "vanillakernel" && nobs > np) {
    ## linear DWD
      ftran = "ldwd"
      if(abs(qval - as.integer(qval)) < .Machine$double.eps) {
        ftran = paste0(ftran, "int")
        qval = as.integer(qval)
      }
      fit = .Fortran(ftran, qval, as.double(x), 
        nobs, np, as.double(y), nlam, ulam, eps, maxit, gam,
        anlam=integer(1), npass=integer(nlam), jerr=integer(1), 
        btmat=double((np + 1) * nlam), PACKAGE="kerndwd")
      anlam = fit$anlam
      alpha = matrix(fit$btmat[seq((np + 1) * anlam)], 
        np + 1, anlam) 
    } else {
    ## kernel DWD
      ftran = "kdwd"
      Kmat = kernelMatrix(kern, x)
      if(abs(qval - as.integer(qval)) < .Machine$double.eps) {
        ftran = paste0(ftran, "int")
        qval = as.integer(qval)
      }
      fit = .Fortran(ftran, qval, as.double(Kmat), 
        nobs, as.double(y), nlam, ulam, eps, maxit, gam,
        anlam=integer(1), npass=integer(nlam), jerr=integer(1), 
        alpmat=double((nobs + 1) * nlam), PACKAGE="kerndwd")
      anlam = fit$anlam
      alpha = matrix(fit$alpmat[seq((nobs + 1) * anlam)], 
        nobs + 1, anlam) 
    }
  } else {
    if (length(wt) != nobs)
      stop("The length of the weight vector is not n.")
    if (any(wt < 0))
      stop("The weights must be nonnegative.")
    if (class(kern)[[1]] == "vanillakernel" && nobs > np) {
    ## weighted linear DWD
      ftran = "wldwd"
      if(abs(qval - as.integer(qval)) < .Machine$double.eps) {
        ftran = paste0(ftran, "int")
        qval = as.integer(qval)
      }
      fit = .Fortran(ftran, qval, as.double(x), as.double(wt),
        nobs, np, as.double(y), nlam, ulam, eps, maxit, gam,
        anlam=integer(1), npass=integer(nlam), jerr=integer(1), 
        btmat=double((np + 1) * nlam), PACKAGE="kerndwd")
      anlam = fit$anlam
      alpha = matrix(fit$btmat[seq((np + 1) * anlam)], 
        np + 1, anlam) 
    } else {
    ## weighted kernel DWD
      ftran = "wkdwd"
      Kmat = kernelMatrix(kern, x)
      if(abs(qval - as.integer(qval)) < .Machine$double.eps) {
        ftran = paste0(ftran, "int")
        qval = as.integer(qval)
      }
      fit = .Fortran(ftran, qval, as.double(Kmat), as.double(wt),
        nobs, as.double(y), nlam, ulam, eps, maxit, gam,
        anlam=integer(1), npass=integer(nlam), jerr=integer(1), 
        alpmat=double((nobs + 1) * nlam), PACKAGE="kerndwd")
      anlam = fit$anlam
      alpha = matrix(fit$alpmat[seq((nobs + 1) * anlam)], 
        nobs + 1, anlam)  
    } 
  }
  ####################################################################
  ## wrap up output
  info = list(qval = qval, eps = eps, maxit = signif(maxit),
    kern = capture.output(show(kern)))
  if (!missing(wt)) info = c(info, wt=list(wt))
  outlist = list(alpha = alpha, lambda = ulam[seq(anlam)], 
    npass = fit$npass, jerr = fit$jerr, info = info)
  class(outlist) = c("dwdpath")
  outlist
}

lumpath = function(x, y, nobs, np, kern, aval, cval, 
  ulam, nlam, wt, eps, maxit, gam) {
  #################################################################### 
  ## check index 
  if (aval <= 0 || cval <=0) stop("aval and cval must be positive.")
  aval = as.double(aval)
  cval = as.double(cval)
  ####################################################################
  if (missing(wt) || is.null(wt)) {
    if (class(kern)[[1]] == "vanillakernel" && nobs > np) {
    ## linear LUM
      ftran = "llum"
      if(abs(aval - as.integer(aval)) < .Machine$double.eps) {
        ftran = paste0(ftran, "int")
        aval = as.integer(aval)
      }
      fit = .Fortran(ftran, aval, cval, as.double(x), 
        nobs, np, as.double(y), nlam, ulam, eps, maxit, gam,
        anlam=integer(1), npass=integer(nlam), jerr=integer(1), 
        btmat=double((np + 1) * nlam), PACKAGE="kerndwd")
      anlam = fit$anlam
      alpha = matrix(fit$btmat[seq((np + 1) * anlam)], 
        np + 1, anlam) 
    } else {
    ## kernel LUM
      ftran = "klum"
      if(abs(aval - as.integer(aval)) < .Machine$double.eps) {
        ftran = paste0(ftran, "int")
        aval = as.integer(aval)
      }
      Kmat = kernelMatrix(kern, x)
      fit = .Fortran(ftran, aval, cval, as.double(Kmat), 
        nobs, as.double(y), nlam, ulam, eps, maxit, gam,
        anlam=integer(1), npass=integer(nlam), jerr=integer(1), 
        alpmat=double((nobs + 1) * nlam), PACKAGE="kerndwd")
      anlam = fit$anlam
      alpha = matrix(fit$alpmat[seq((nobs + 1) * anlam)], 
        nobs + 1, anlam) 
    }
  } else {
    if (length(wt) != nobs)
      stop("The length of the weight vector is not n.")
    if (any(wt < 0))
      stop("The weights must be nonnegative.")
    if (class(kern)[[1]] == "vanillakernel" && nobs > np) {
    ## weighted linear LUM
      ftran = "wllum"
      if(abs(aval - as.integer(aval)) < .Machine$double.eps) {
        ftran = paste0(ftran, "int")
        aval = as.integer(aval)
      }
      fit = .Fortran(ftran, aval, cval, as.double(x), as.double(wt),
        nobs, np, as.double(y), nlam, ulam, eps, maxit, gam,
        anlam=integer(1), npass=integer(nlam), jerr=integer(1), 
        btmat=double((np + 1) * nlam), PACKAGE="kerndwd")
      anlam = fit$anlam
      alpha = matrix(fit$btmat[seq((np + 1) * anlam)], 
        np + 1, anlam) 
    } else {
    ## weighted kernel LUM
      ftran = "wklum"
      if(abs(aval - as.integer(aval)) < .Machine$double.eps) {
        ftran = paste0(ftran, "int")
        aval = as.integer(aval)
      }
      Kmat = kernelMatrix(kern, x)
      fit = .Fortran(ftran, aval, cval, as.double(Kmat), as.double(wt),
        nobs, as.double(y), nlam, ulam, eps, maxit, gam,
        anlam=integer(1), npass=integer(nlam), jerr=integer(1), 
        alpmat=double((nobs + 1) * nlam), PACKAGE="kerndwd")
      anlam = fit$anlam
      alpha = matrix(fit$alpmat[seq((nobs + 1) * anlam)], 
        nobs + 1, anlam)  
    } 
  }
  ####################################################################
  ## wrap up output
  info = list(aval = aval, cval=cval, eps = eps, 
    maxit = signif(maxit), kern = capture.output(show(kern)))
  if (!missing(wt)) info = c(info, wt=list(wt))
  outlist = list(alpha = alpha, lambda = ulam[seq(anlam)], 
    npass = fit$npass, jerr = fit$jerr, info = info)
  class(outlist) = c("lumpath")
  outlist
}

otherpath = function(x, y, nobs, np, kern, const, dev,
  ulam, nlam, wt, eps, maxit, gam) {
  ####################################################################
  if (missing(const)) 
    stop("When 'family = other', 'const' has to be specified.")
  dev = match.fun(dev)
  if (missing(wt) || is.null(wt)) wt = rep(1, nobs) else {
    if (length(wt) != nobs)
      stop("The length of the weight vector is not n.")
    if (any(wt < 0))
      stop("The weights must be nonnegative.")
  }
  if (class(kern)[[1]] == "vanillakernel" && nobs > np)
    fit = otherlinrpath(x, y, nobs, np, kern, const, dev,
      ulam, nlam, wt, eps, maxit, gam) else
    fit = otherkernpath(x, y, nobs, np, kern, const, dev,
      ulam, nlam, wt, eps, maxit, gam)
  class(fit) = c("otherpath")
  fit
}

otherkernpath = function(x, y, nobs, np, kern, const, dev,
      ulam, nlam, wt, eps, maxit, gam) {
  Kmat = kernelMatrix(kern, x)
  WK = sweep(Kmat, MARGIN=1, wt, '*')
  npass = rep(0, nlam)
  jerr = 0
  anlam = 0
  alp = rep(0, nobs + 1)
  alpmat = matrix(NA, nobs + 1, nlam)
  Plam = matrix(NA, (nobs + 1), (nobs + 1))
  Plam[1, 1] = sum(wt)
  WKs = colSums(WK)
  Plam[1, -1] = WKs
  Plam[-1, 1] = WKs
  Plam[-1, -1] = Kmat %*% WK
  Plam = Plam + gam * diag(nobs + 1)
  Plaml = Plam
  minv = 1/const
  rvec = zvec = ds = rep(0, nobs)
  for (l in seq.int(nlam)) {
    Plaml[-1, -1] = Plam[-1, -1] + 2 * nobs * ulam[l] * minv * Kmat
    Pinv = solve(Plaml)
    while (1) {
      for (i in seq.int(nobs)) ds[i] = dev(rvec[i])
      zvec = wt * y * ds
      ru = sum(zvec)
      rd = Kmat %*% (zvec + 2 * nobs * ulam[l] * alp[-1])
      dif = - minv * Pinv %*% c(ru, rd)
      alp = alp + dif
      rvec = rvec + y * (dif[1] + Kmat %*% dif[-1])
      npass[l] = npass[l] + 1
      if (sum(dif * dif) < eps) break
      if (sum(npass) > maxit) break
    }
    if (sum(npass) > maxit) {
      jerr = -l
      break
    }
    anlam = l
    alpmat[, l] = alp
  }
  alpha = alpmat[, seq(anlam)]
  ####################################################################
  ## wrap up output
  info = list(eps = eps, dev = dev,
    maxit = signif(maxit), kern = capture.output(show(kern)), 
    wt = wt)
  outlist = list(alpha = alpha, lambda = ulam[seq(anlam)], 
    npass = npass, jerr = jerr, info = info)
  outlist
}

otherlinrpath = function(x, y, nobs, np, kern, const, dev,
      ulam, nlam, wt, eps, maxit, gam) {
  bt = rep(0, np + 1)
  alpmat = matrix(NA, np + 1, nlam)
  minv = 1 / const
  npass = rep(0, nlam)
  jerr = 0
  anlam = 0
  WX = sweep(x, MARGIN=1, wt, '*')
  WXsum = colSums(WX)
  Plam = matrix(NA, (np + 1), (np + 1))
  Plam[1, 1] = sum(wt)
  Plam[1, -1] = WXsum
  Plam[-1, 1] = WXsum
  Plam[-1, -1] = crossprod(x, WX)
  Plam = Plam + gam * diag(np + 1)
  Plaml = Plam
  rvec = zvec = ds = rep(0, nobs)
  for (l in seq.int(nlam)) {
    diag(Plaml)[-1] = diag(Plam)[-1] + 2 * nobs * ulam[l] * minv
    Pinv = solve(Plaml)
    rd = rep(NA, np + 1)
    while (1) {
      for (i in seq.int(nobs)) ds[i] = dev(rvec[i])
      zvec = wt * y * ds
      rd[1] = sum(zvec)
      rd[-1] = 2 * nobs * ulam[l] * bt[-1] + crossprod(x, zvec)
      dif = - minv * Pinv %*% rd
      bt = bt + dif
      rvec = rvec + y * (dif[1] + x %*% dif[-1])
      npass[l] = npass[l] + 1
      if (sum(dif * dif) < eps) break
      if (sum(npass) > maxit) break
    }
    if (sum(npass) > maxit) {
      jerr = -l
      break
    }
    anlam = l
    alpmat[, l] = bt
  }
  alpha = alpmat[, seq(anlam)]
  ####################################################################
  ## wrap up output
  info = list(eps = eps, dev = dev,
    maxit = signif(maxit), kern = capture.output(show(kern)), 
    wt = wt)
  outlist = list(alpha = alpha, alpha = alpha, lambda = ulam[seq(anlam)], 
    npass = npass, jerr = jerr, info = info)
  outlist
}


