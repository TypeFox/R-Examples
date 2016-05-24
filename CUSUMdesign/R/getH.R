getH = function(distr=NULL, ARL=NULL, ICmean=NULL, ICsd=NULL, 
    OOCmean=NULL, OOCsd=NULL, ICprob=NULL, OOCprob=NULL, ICvar=NULL, 
    IClambda=NULL, samp.size=NULL, ref=NULL, winsrl=NULL, winsru=NULL, 
    type=c("fast initial response", "zero start", "steady state")) {
  if (is.null(distr)) {
    cat('1 = Normal location, \n 
         2 = Normal variance, \n
         3 = Poisson, \n 
         4 = Binomial, \n 
         5 = Negative binomial, \n
         6 = Inv Gaussian mean. \n')
    stop("Sepcify a distribution.")
  }
  distrs = c('Normal location','Normal variance',
    'Poisson', 'Binomial', 'Negative binomial',
    'Inv Gaussian mean')
  nonmen = length(distrs)
  if (!(distr %in% seq(nonmen))) stop ("'distr' is not a valid choice.")
  
  ## parameter settings
  eps0 = 0.05
  eps = 0.001
  int10 = as.integer(10)
  plus = 1.0
  denrat = 0.0
  isint = FALSE
  step0 = 0.25
  quar3 = 0.75
  quar = 0.25
  arg2 = 0.0
  intone = as.integer(1)

  ## Winsorizing constants
  if (is.null(winsrl)) winsrl = -999 else winsrl = as.double(winsrl)
  if (is.null(winsru)) winsru = 999 else winsru = as.double(winsru)
  ## zero-start (Z), steady state (S), or fast initial response (F)
  if (is.null(type)) {
    cat ("type is missing. Set type = 'F' (FIR). \n")
    type = "F"
  }
  type = tolower(type)
  type = match.arg(type)
  runtype = match(type, 
    c("zero start", "fast initial response", "steady state"))
  if (is.null(ARL)) stop ("ARL is missing.") else ARL = as.double(ARL)
  
  ## read parameters for different distributions
  if (distr == 1) {
    # Normal location
    if (is.null(ICmean)) stop ("In-control mean 'ICmean' is missing.")
    if (is.null(ICsd)) stop ("In-control sd 'ICsd' is missing.")
    if (is.null(OOCmean)) stop ("Out-of-control mean 
      'OOCmean' is missing.")
    amu = as.double(ICmean)
    amu1 = as.double(OOCmean)
    sigma = as.double(ICsd)
    if (is.null(ref)) 
      ref = 0.5 * abs(amu1 - amu) / sigma else
      cat("ref is not user-specifed for normal location, \n")
    offset = 0.5 * (amu + amu1)
    cat(sprintf("The reference value is %12.3f \n", offset))
    ndis = intone
  }

  if (distr == 2) {
    ## Normal variance
    if (is.null(samp.size)) stop ("Rational group size 'samp.size' 
      is missing.") else nsamp = as.integer(samp.size)
    if (nsamp < 2) {
      cat("Using squared deviations from true mean. \n")
      nsamp = as.integer(2)
    }
    integ = nsamp - intone
    # change the common variable 'ndf' in varup & vardn subroutines
    .Fortran("cmpar2", integ)
    if (is.null(ICsd)) stop ("In-control sd 'ICsd' is missing.")
    if (is.null(OOCsd)) stop ("Out-of-control sd 'OOCsd' is missing.")
    sigma0 = as.double(ICsd)
    sigma1 = as.double(OOCsd)
    varrat = (sigma0 / sigma1) ^ 2
    if (is.null(ref))
      ref = log(varrat) / (varrat - 1) else
      cat("ref is not user-specifed for normal variance, \n")
    offset = ref * sigma0 * sigma0
    cat(sprintf("The reference value is %12.3f \n", offset))
    ndis = as.integer(2)
    if (sigma1 < sigma0) {
      plus = -1
      ndis = as.integer(3)
    }
  }
  
  if (distr == 3) {
    ## Poisson
    if (is.null(ICmean)) stop ("In-control mean 'ICmean' is missing.")
    if (is.null(OOCmean)) stop ("Out-of-control mean 
      'OOCmean' is missing.")
    arg2 = as.double(ICmean)
    amu1 = as.double(OOCmean)
    if (min(arg2, amu1) <= 0) 
      stop("Invalid - means must be strictly positive.'")
    if (is.null(ref))
      ref = (amu1 - arg2) / log(amu1 / arg2) else
      ref = as.double(ref)
    cat(sprintf("The reference value is %12.3f \n", ref))
    ndis = as.integer(4)
    if (amu1 < arg2) {
      plus = -1
      ndis = as.integer(5)
    }
    isint = TRUE
    step0 = as.integer(arg2 * 0.5 + 1) 
  }

  if (distr == 4) {
    ## Binomial
    if (is.null(samp.size)) stop ("Sample size 'samp.size' is missing.")
    if (is.null(ICprob)) stop ("In-control probability
      'ICprob' is missing.")
    if (is.null(OOCprob)) stop ("Out-of-control probability
      'OOCprob' is missing.")
    integ = as.integer(samp.size)
    # change the common variable 'nbig' in binup & bindn subroutines
    .Fortran("cmpar2", integ)
    arg2 = as.double(ICprob)
    pi1 = as.double(OOCprob)
    if(min(pi1, arg2) <= 0 || max(pi1, arg2) > 1)
      stop("Invalid IC or OOC probability.")
    if (is.null(ref))
      ref = -integ * log((1 - pi1)/(1 - arg2)) / log(pi1 * (1 - arg2)/
      (arg2 * ((1 - pi1)))) else
      ref = as.double(ref)
    cat(sprintf("The reference value is %12.3f \n", ref))
    ndis = as.integer(6)
    if (pi1 < arg2) {
      plus = -1
      ndis = as.integer(7)
    } 
    isint = TRUE
    step0 = as.integer((integ * arg2) * 0.5 + 1)
  }

  if (distr == 5) {
    ## Negative binomial
    if (is.null(ICmean)) stop ("In-control mean 'ICmean' is missing.")
    if (is.null(ICvar)) stop ("In-control variance 'ICvar' is missing.")
    if (is.null(OOCmean)) stop ("Out-of-control mean
      'OOCmean' is missing.")
    aver = as.double(ICmean)
    varian = as.double(ICvar)
    aver1 = as.double(OOCmean)
    if(min(aver, aver1, varian - aver) <= 0) 
      stop("Means must be > 0 and variance > in-control mean.")
    arg2 = aver / (varian - aver)
    realno = aver * arg2
    # change the common variable 'r' in nbinup & nbindn subroutines
    .Fortran("cmpar", realno)
    c1 = arg2 * aver / aver1
    if (is.null(ref))
      ref = - realno * log((c1*(1+arg2))/(arg2*(1+c1))) /
        log((1+arg2)/(1+c1)) else
      ref = as.double(ref)
    cat(sprintf("The reference value is %12.3f \n", ref))
    ndis = as.integer(8)
    if (aver1 < aver) {
      plus = -1
      ndis = as.integer(9)
    }
    isint = TRUE
    step0 = as.integer(aver * 0.5 + 1)
  }

  if (distr == 6) {
  ## Inverse Gaussian mean
    if (is.null(ICmean)) stop ("In-control mean 'ICmean' is missing.")
    if (is.null(IClambda)) stop ("In-control lambda 
      'IClambda' is missing.")
    if (is.null(OOCmean)) stop ("Out-of-control mean 'OOCmean' is missing.")
    arg2 = as.double(ICmean)
    realno = as.double(IClambda)
    # change the common variable 'alam' in gauiup & gauidn subroutines
    .Fortran("cmpar", realno)
    amu1 = as.double(OOCmean)
    if (min(arg2, realno, amu1) <= 0.0)
      stop("All parameters must be strictly positive.")
    if (is.null(ref))
      ref = 2 * arg2 * amu1 / (arg2 + amu1) else
      ref = as.double(ref)
    cat(sprintf("The reference value is %12.3f \n", ref))
    step0 = min(arg2, realno, amu1) * 0.5
    ndis = as.integer(10)
    if (amu1 < arg2){
      plus = -1
      ndis = as.integer(11)
    }
  }

  ## Do preliminary ranging
  maxfau = as.integer(0)
  hlo = 0
  ihlo = as.integer(0)
  ihhi = as.integer(0)
  intval = as.integer(0)
  istep = as.integer(0)
  alow = 1
  step = step0
  temp = 1

  ## If integer - find rational denominator
  if (isint) {
    den = .Fortran('getden', ref=as.double(ref), temp=as.double(temp), 
      denrat=as.double(denrat), as.logical(FALSE), as.logical(FALSE))
    ref = den$ref
    temp = den$temp
    denrat = den$denrat
    istep = step * denrat
    istep = min(istep, int10)
  }
  
  ## Ranging
  for (inner in seq(40)) {
    step = step * 1.2
    hhi = hlo + step
    if (isint) {
      istep = min(int10, istep + 2)
      ihhi = ihlo + istep
      if (ihhi > 256) 
        stop("Can not handle that K. Try a simpler one")
      hhi = as.double(ihhi) / denrat
    }
    winlo = winsrl
    winhi = winsru
    if (plus < 0) {
      winlo = -winsru
      winhi = -winsrl
    }
    cal = .Fortran('calcus', as.double(hhi), as.double(plus*ref), 
      as.integer(ndis), as.double(denrat), as.double(arg2), 
      as.double(winlo), as.double(winhi), regarl=as.double(0), 
      firarl=as.double(0), ssarl=as.double(0), eps0=eps0, 
      esterr=as.double(0), ifault=as.integer(0))
    regarl = cal$regarl
    firarl = cal$firarl
    ssarl = cal$ssarl
    esterr = cal$esterr
    ifault = cal$ifault
    maxfau = max(maxfau,ifault)
    if (ifault != 0) {
      good = 0
      if (esterr > 0) good = -log10(esterr) 
      cat(sprintf(' h %11.4f arls %9.1f %7.1f %7.1f \n', 
        hhi, regarl, firarl, ssarl))
    }
    gotarl = switch(runtype, regarl, firarl, ssarl)
    ahi = gotarl
    if (gotarl > ARL) break
    hlo = hhi
    ihlo = ihhi
    alow = gotarl
  }
  if (gotarl <= ARL) stop ("Ranging failed. Exit.")

  ## Refinement
  logscl = (min(alow, ahi, ARL) > 0)
  if (logscl) {
    alowlg = log(alow)
    ahilg = log(ahi)
    arllg = log(ARL)
  } else {
    alowlg = alow
    ahilg = ahi
    arllg = ARL 
  }

  for (inner in seq(20)) {
    ratio = (arllg - alowlg) / (ahilg - alowlg)
    ratio = min(quar3, max(quar, ratio))
    test = hlo + ratio * (hhi - hlo)
    if (isint && ((ihhi - ihlo) == 1)) {
      cat(sprintf("k %8.4f   h %8.4f  ARL %10.2f  
             h %8.4f  ARL %10.2f \n", ref, hlo, alow, hhi, ahi))
      test = hhi
      gothi = ahi
      break
    }
    if (isint) {
      intval = as.integer((ihlo + ihhi + intone) / 2)
      intval = max(ihlo + intone, 
        min(ihhi - intone, intval))
      test = as.double(intval) / denrat
    } 
    cal2 = .Fortran('calcus', as.double(test), as.double(plus*ref), 
      as.integer(ndis), as.double(denrat), as.double(arg2), 
      as.double(winlo), as.double(winhi), regarl=as.double(0), 
      firarl=as.double(0), ssarl=as.double(0), eps=eps, 
      esterr=as.double(1), ifault=as.integer(0))
    regarl = cal2$regarl
    firarl = cal2$firarl
    ssarl = cal2$ssarl
    esterr = cal2$esterr
    ifault = cal2$ifault
    maxfau = max(maxfau,ifault)
    if (ifault == 0) {
#      cat(sprintf(' h %11.4f arls %9.1f %7.1f %7.1f \n', 
#        test, regarl, firarl, ssarl))
    } else {
      good = -log10(esterr)
      cat(sprintf(' h %11.4f arls %9.1f %7.1f %7.1f %7.1f 
        estimated good digits \n', 
        test, regarl, firarl, ssarl, good))
    }
    gotarl = switch(runtype, regarl, firarl, ssarl)
    if (gotarl < ARL) {
      hlo = test
      alowlg = gotarl
      if (logscl) alowlg = log(gotarl)
      ihlo = intval
      alow = gotarl
    } else {
      hhi = test
      ahilg = gotarl
      if (logscl) ahilg = log(gotarl)
      ihhi = intval
      ahi = gotarl
    }
    if (abs(gotarl / ARL - 1.0) < eps) break
  }
  if (inner > 40) {
    stop(sprintf("Convergence failed. Best guesses.
      \n k %8.4f   h %8.4f ARL %10.2f 
        h %8.4f ARL %10.2f \n", ref, hlo, alow, hhi, ahi))
  }
  if(ndis == 1) {
    ## Normal location
    fit = .Fortran('calcus', as.double(test), as.double(-abs(ref)), 
      as.integer(ndis), as.double(denrat), as.double(arg2), 
      as.double(winlo), as.double(winhi),
      regarl=as.double(0), firarl=as.double(0), ssarl=as.double(0),
      eps, esterr=as.double(1), ifault=as.integer(0))
    test = test * sigma
    cat(sprintf("DI %11.3f IC ARL %9.1f  OOC ARL Zero start %7.1f
      FIR %7.1f SS %7.1f \n", test, gotarl, fit$regarl,
      fit$firarl, fit$ssarl)) 
    res = list(DI = test, ref = offset, 
      IC_ARL = gotarl, OOCARL_Z = fit$regarl, 
      OOCARL_F = fit$firarl, OOCARL_S = fit$ssarl)
  }

  if(ndis == 2 || ndis == 3) {
    ## Normal variance
    fit = .Fortran('calcus', as.double(test*varrat), 
      as.double(plus * ref * varrat), 
      as.integer(ndis), as.double(denrat), as.double(arg2), 
      as.double(winlo), as.double(winhi),
      regarl=as.double(0), firarl=as.double(0), ssarl=as.double(0),
      eps, esterr=as.double(1), ifault=as.integer(0))
    test = test * sigma0 * sigma0
    cat(sprintf("DI %11.3f IC ARL %9.1f  OOC ARL Zero start %7.1f
      FIR %7.1f SS %7.1f \n", test, gotarl, fit$regarl,
      fit$firarl, fit$ssarl)) 
    res = list(DI = test, ref = offset, 
      IC_ARL = gotarl, OOCARL_Z = fit$regarl, 
      OOCARL_F = fit$firarl, OOCARL_S = fit$ssarl)
  }

  if(ndis == 4 || ndis == 5) {
    ## Poisson
    fit = .Fortran('calcus', as.double(test), 
      as.double(plus * ref), 
      as.integer(ndis), as.double(denrat), as.double(amu1), 
      as.double(winlo), as.double(winhi),
      regarl=as.double(0), firarl=as.double(0), ssarl=as.double(0),
      eps, esterr=as.double(1), ifault=as.integer(0))
    cat(sprintf("DI %11.3f IC ARL %9.1f  OOC ARL Zero start %7.1f
      FIR %7.1f SS %7.1f \n", test, gothi, fit$regarl,
      fit$firarl, fit$ssarl)) 
    res = list(DI = test, ref = ref, 
      IC_ARL = gotarl, OOCARL_Z = fit$regarl, 
      OOCARL_F = fit$firarl, OOCARL_S = fit$ssarl)
  }

  if(ndis == 6 || ndis == 7) {
    ## Bionomial
    fit = .Fortran('calcus', as.double(test), 
      as.double(plus * ref), 
      as.integer(ndis), as.double(denrat), as.double(pi1), 
      as.double(winlo), as.double(winhi),
      regarl=as.double(0), firarl=as.double(0), ssarl=as.double(0),
      eps, esterr=as.double(1), ifault=as.integer(0))
    cat(sprintf("DI %11.3f IC ARL %9.1f  OOC ARL Zero start %7.1f
      FIR %7.1f SS %7.1f \n", test, gothi, fit$regarl,
      fit$firarl, fit$ssarl)) 
    res = list(DI = test, ref = ref, 
      IC_ARL = gotarl, OOCARL_Z = fit$regarl, 
      OOCARL_F = fit$firarl, OOCARL_S = fit$ssarl)
  }

  if(ndis == 8 || ndis == 9) {
    ## Negative Bionomial
    fit = .Fortran('calcus', as.double(test), 
      as.double(plus * ref), 
      as.integer(ndis), as.double(denrat), as.double(c1), 
      as.double(winlo), as.double(winhi),
      regarl=as.double(0), firarl=as.double(0), ssarl=as.double(0),
      eps, esterr=as.double(1), ifault=as.integer(0))
    cat(sprintf("DI %11.3f IC ARL %9.1f  OOC ARL Zero start %7.1f
      FIR %7.1f SS %7.1f \n", test, gothi, fit$regarl,
      fit$firarl, fit$ssarl)) 
    res = list(DI = test, ref = ref, 
      IC_ARL = gotarl, OOCARL_Z = fit$regarl, 
      OOCARL_F = fit$firarl, OOCARL_S = fit$ssarl)
  }

   if(ndis == 10 || ndis == 11) {
    ## Inverse Gaussian mean
    fit = .Fortran('calcus', as.double(test), 
      as.double(plus * ref), 
      as.integer(ndis), as.double(denrat), as.double(amu1), 
      as.double(winlo), as.double(winhi),
      regarl=as.double(0), firarl=as.double(0), ssarl=as.double(0),
      eps, esterr=as.double(1), ifault=as.integer(0))
    cat(sprintf("DI %11.3f IC ARL %9.1f  OOC ARL Zero start %7.1f
      FIR %7.1f SS %7.1f \n", test, gotarl, fit$regarl,
      fit$firarl, fit$ssarl)) 
    res = list(DI = test, ref = ref, 
      IC_ARL = gotarl, OOCARL_Z = fit$regarl, 
      OOCARL_F = fit$firarl, OOCARL_S = fit$ssarl)
  }
  ## Output the result
  res
}


