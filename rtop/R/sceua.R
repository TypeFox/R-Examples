


sceua = function(OFUN, pars, lower, upper, maxn = 10000, kstop = 5, pcento = 0.01,
    ngs = 5, npg = 5, nps = 5, nspl = 5, mings = 5, iniflg = 1, iprint = 0, iround = 3, 
    peps = 0.0001, plog = rep(FALSE,length(pars)), implicit = NULL, ...) {

# OFUN - objective function
# pars - starting values
# lower - lower bounds
# upper - upper bounds
# maxn - maximum number of iterations
# kstop - number of shuffling loops in which the criterion value must change 
#         by the given percentage before optimization is terminated
# pcento - percentage by which the criterion value must change in given number of shuffling loops
# ngs - number of complexes in the initial population
# npg = number of points in each complex
# npt = total number of points in initial population (npt=ngs*npg)
# nps = number of points in a sub-complex
# nspl = number of evolution steps allowed for each complex before
#     complex shuffling
# mings = minimum number of complexes required, if the number of
#     complexes is allowed to reduce as the optimization proceeds
# iniflg = flag on whether to include the initial point in population
#     = 0, not included
#     = 1, included
# iprint = flag for controlling print-out after each shuffling loop
#     = 0, print information on the best point of the population
#      = 1, print information on every point of the population
# implicit = function for implicit boundaries (e.g. sum(par[4]+par[5]) < 1)

  npars = length(pars)
  if (length(plog) == 1) 
    plog = rep(plog,npars) 
  if (length(upper) != npars | length(lower) != npars | length(plog) != npars) 
     stop("pars, upper, lower and plog must be of same length, plog can alternatively be of length 1")
  pars = ifelse(plog,log10(pars),pars)
  upper = ifelse(plog,log10(upper),upper)
  lower = ifelse(plog,log10(lower),lower)

  nloop = 0
  npt = ngs * npg
  loop = 0
  bound = upper-lower
  criter = rep(1e10,20)
  parset = matrix(nrow = npt,ncol = npars)
  xf = rep(1e10,npt)
  icall = 1

  lpars = ifelse(plog,10^pars,pars)
  fa = OFUN(lpars,...)
  if (iprint > 0 && icall %% iprint == 0) cat(icall,signif(fa,iround), "\n")
  parset[1,] = pars
  xf[1] = fa
  stdinit = rep(1,npars)
  for (ii in ifelse(iniflg == 1,2,1):npt) {
    parset[ii,] = getpnt(idist = 1,lower,upper,stdinit,lower, implicit)
    lpars = ifelse(plog,10^parset[ii,],parset[ii,])
    xf[ii] = OFUN(lpars,...)
    icall = icall + 1
    if (iprint > 0 && icall %% iprint == 0) cat(icall,round(xf[ii],iround), "\n")
  }

  parset = parset[order(xf),]
  xf = sort(xf)
  bestpar = parset[1,]
  worstpar = parset[npt,]
  bestf = xf[1]
  worstf = xf[npt]

  parsttout = parstt(npt,npars,parset,bound, peps)
  ipcnvg = parsttout$ipcnvg
  gnrng = parsttout$gnrng
  parstd = parsttout$parstd
  
  while (TRUE) {
    nloop = nloop + 1
    for (igs in 1:ngs) {
      karr = (c(1:npg)-1)*ngs + igs
      cx = parset[karr,]
      cf = xf[karr]
      for (loop in 1:nspl) {
        kpos = 1
        if (nps == npg) {
          lcs = c(1:nps)
        } else while(TRUE) {
          lpos = 1+round(npg+.5-sqrt((npg+.5)^2-npg*(npg+1)*runif(1)))
          lcs[kpos] = lpos
          if (sum(duplicated(lcs))==0) kpos = kpos + 1
          if (kpos > nps) break 
        }
        lcs = sort(lcs)
        s = cx[lcs,]
        sf = cf[lcs]  
        cceout = cce(OFUN,npars,nps,s,sf,lower,upper,parstd,icall,maxn,iprint,iround,bestf,plog,implicit,...)
        s = cceout$s
        sf = cceout$sf
        icall = cceout$icall
        cx[lcs,] = s
        cf[lcs] = sf
        cx = cx[order(cf),]
        cf = sort(cf)
      }
      parset[karr,] = cx
      xf[karr] = cf
    }
    parset = parset[order(xf),]
    xf = sort(xf)
    bestpar = parset[1,]
    worstpar = parset[npt,]
    bestf = xf[1]
    worstf = xf[npt]
    parsttout = parstt(npt,npars,parset,bound, peps)
    ipcnvg = parsttout$ipcnvg
    gnrng = parsttout$gnrng
    parstd = parsttout$parstd
    fbestf = criter[kstop]
    concrit = 2*(fbestf-bestf)/(fbestf+bestf)
    criter[2:length(criter)] = criter[1:(length(criter)-1)]
    criter[1] = bestf
    if (iprint >= 0) cat(icall,"best",
        signif(bestf,iround), "function convergence", signif(concrit,iround)/pcento,
        "parameter convergence", gnrng/peps, "\n")

    if (concrit < pcento & ipcnvg == 1) break
    if (icall > maxn) break
    if (ngs > mings) {
      compout = comp(npars,npt,ngs,npg,parset,xf)
      ngs = ngs -1
      parset = compout$parset
      xf = compout$xf
    }
  }
  bestpar = ifelse(plog,10^bestpar,bestpar)
  return(list(par = bestpar, value = xf[1], convergence = list(funConvergence = signif(concrit,iround)/pcento, parConvergence = gnrng/peps),
      counts = icall, iterations = nloop))
}

comp = function(npars,npt,ngs,npg,parset,xf){
  xn = parset
  xfn = xf
  for (igs in 1:ngs) {
    karr1 = (c(1:npg)-1)*ngs + igs
    karr2 = (c(1:npg)-1)*(ngs-1) + igs
    xn[karr2,] = parset[karr1,]
    xfn[karr2] = xf[karr1]
  }
  return(list(parset = xn,xf = xfn))
}

cce = function(OFUN,npars,nps,s,sf,lower,upper,parstd,icall,maxn,iprint,iround,bestf,plog, implicit,...) {
  alpha = 1.
  beta = 0.5
  n = dim(s)[1]
  sb = s[1,]
  sw = s[n,]
  ce = colMeans(s)
  fw = sf[n]
  snew = ce+alpha*(ce-sw)
#  print(icall)
  if (chkcst(snew,lower,upper, implicit) >0) snew = getpnt(2, lower, upper, parstd, sb, implicit)
#  print(snew)
  lpars = ifelse(plog,10^snew,snew)
  fnew = OFUN(lpars,...)
  icall = icall + 1
  if (iprint > 0 && icall %% iprint == 0) cat(icall, signif(fnew, iround), signif(bestf, iround), "\n")
  if (fnew > fw) {
    snew = ce-beta*(ce-sw)
    lpars = ifelse(plog,10^snew,snew)
    fnew = OFUN(lpars,...)
    icall = icall + 1
    if (iprint > 0 && icall %% iprint == 0) cat(icall, signif(fnew, iround), signif(bestf, iround), "\n")
    if (fnew > fw) {
      snew = getpnt(2,lower,upper,parstd,sb, implicit)
      lpars = ifelse(plog,10^snew,snew)
      fnew = OFUN(lpars,...)
      icall = icall + 1
      if (iprint > 0 && icall %% iprint == 0) cat(icall, signif(fnew, iround), signif(bestf, iround), "\n")
    }
  }
  s[n,] = snew
  sf[n] = fnew
  return(list(s = s, sf = sf, icall = icall))
}

chkcst = function(parlocal,lower,upper, implicit) {                                                  
 ibound = if (sum(mapply(FUN = function(x,y,z)
               max(y-x,x-z,0),parlocal,lower,upper)) > 0) 1 else 0
 if (ibound == 0 & length(parlocal) >1 & !is.null(implicit)) {
# Possibility to include implicit constraints 
   if (!is.function(implicit)) stop("implicit has to be a function")
   ibound = implicit(parlocal)
 }
 ibound
}


getpnt = function(idist,lower,upper,std,pari, implicit){
#  rand = (ifelse(rep(idist,npars) == 1,runif(npars),rnorm(npars)))
#  print(xi)
#  print(rand)
  ic = 0
  while (TRUE) {
    parj = mapply(FUN = get1p, pari, std = std, lower = lower, upper = upper, 
        MoreArgs = list(idist = idist, implicit = implicit))
    if (chkcst(parj,lower,upper, implicit) == 0) break
    ic = ic + 1
    if (ic > 100) stop("Cannot find a parameter set respecting the fixed or implicit boundaries after 100 iterations")
  }
  return(parj)
}

get1p = function(pari,std,lower,upper,idist, implicit) {
#  print(paste(xi,std,rand,lower,upper))
  while (TRUE) {
    rand = ifelse(idist == 1,runif(1),rnorm(1))
    parj = pari+std*rand*(upper-lower)
#    print(x)
#    print(chkcst(x,lower,upper))
    if (chkcst(parj,lower,upper, implicit) == 0) break
#    print(acdf)
  }
  return(parj)
}

parstt = function(npt,npars,parset,bound, peps) {
  parstd = apply(parset,MARGIN=2,FUN = function(x) sd(x))/bound
  parmin = apply(parset,MARGIN=2,FUN = function(x) min(x))
  parmax = apply(parset,MARGIN=2,FUN = function(x) max(x))
  gsum = sum(log((parmax-parmin)/bound))
  gnrng = exp(gsum/npars)
  ipcnvg = ifelse(gnrng <=peps, 1, 0)
  return(list(ipcnvg = ipcnvg,gnrng = gnrng,parstd = parstd))
}


#FUN = function(pars,target) (pars[1]*pars[2]*pars[3]-target)^2

#p0 = c(1,1,1)
#upper = c(20,20,20)
#lower = c(-20,-20,-20)
#best = sceua(FUN,p0,bl,bu,plog=FALSE,target = 3.75)
# plog is a logical to define if parameters should be logarithmized or not