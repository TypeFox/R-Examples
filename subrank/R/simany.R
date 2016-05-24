simany <- function (sampsize, subsampsizes, sampnum, nbsafe = 5, fun =5, ...) 
{
  DistTypes = colnames(distance(1:5, 5:1))
  lrs = array(dim = c(sampnum, length(subsampsizes), length(DistTypes)))
  lrs2mean = lrs
  lrsmeanmean = array(dim = c(length(subsampsizes), length(DistTypes)))
  scarcities = matrix(ncol = length(subsampsizes), nrow = sampnum)
  indep=is.double(fun)
  if (!indep) { fun <- match.fun(fun) }
  vraicop = list()
  for (s in 1:length(subsampsizes)) {
    if (!indep) {
      simdata = fun(sampsize * sampnum, ...)
      # print(simdata)
      dimension=dim(simdata)[2]
      print(dimension)
      simdata=as.numeric(simdata)
      nboot = nbsafe * subsampsizes[s]^dimension
      vraicoptemp = corc0(simdata, sampsize * sampnum, 
                          dimension, subsampsizes[s], nboot, 42)
      tailcop=subsampsizes[s]^dimension
      nbootreel = vraicoptemp[tailcop + 2]
      vraicop[[s]] = vraicoptemp[1:tailcop]/(nbootreel*subsampsizes[s])
      lrsmeanmean[s,] = distance(vraicop[[s]], rep(1/tailcop, tailcop))
    }
    else {
      dimension=fun
      tailcop=subsampsizes[s]^dimension
      nboot = nbsafe * tailcop
      vraicop[[s]] = rep(1/tailcop, tailcop)
      lrsmeanmean[s,]=0
    }
  }
  for (e in 1:sampnum) {
    for (s in 1:length(subsampsizes)) {
      if (!indep)
        { simdata = fun(sampsize, ...) } else
        { simdata = rnorm(sampsize * dimension) }
      tailcop=subsampsizes[s]^dimension
      nboot = nbsafe * tailcop
      cop = corc0(simdata, sampsize, dimension, subsampsizes[s], 
                  nboot, 42)
      nbootreel = cop[tailcop + 2]
      cop = cop[1:tailcop]/(nbootreel*subsampsizes[s])
      lrs[e, s, ] = distance(cop, rep(1/tailcop, tailcop))
      lrs2mean[e, s, ] = distance(cop, vraicop[[s]])
      scarcities[e, s] = sum(cop == 0)/(nbootreel * subsampsizes[s])
    }
  }
  return(list(lrs = lrs, lrs2mean = lrs2mean, lrsmeanmean=lrsmeanmean, scarcities = scarcities, 
              DistTypes = DistTypes))
}