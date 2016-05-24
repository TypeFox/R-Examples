zoomProfile <- function(fixedlist=NA, extrap=F, locchull=NULL,
                      templateinkg=blackbox.getOption("rosglobal")$par, ##but another $par can be passed if this one is not in locchull
                      precision
) { ## Nb in metric units for metric plot.
  ## optimization in a larger space of an objective function that iteratively devalues being far from a subspace
  namesfixed <- names(fixedlist)
  fixedvec <- unlist(fixedlist)[namesfixed]
  ## fixedlist argument was already logscale, if relevant
  locchullmax <- apply(locchull$vertices, 2, max)
  locchullmin <- apply(locchull$vertices, 2, min)
  localmaxrange <- locchullmax-locchullmin
  rangefixed <- localmaxrange[namesfixed] ## in scale of hull.... which must have been determined by islogscale()
  cost <- 1
  newpar <- fromFONKtoanyspace(templateinkg, colnames(locchull$vertices))
  repeat { ## loop starting from the unconstrained maximum with objective function increasingly peaked on the constraint
    oldpar <- newpar
    objfn <- function(x) { ## x is in hull space, hence fixedvec must be too
      tofKpredict.nohull(x, fixedlist=NA)-cost*mean(((x[namesfixed]-fixedvec)/rangefixed)^2)
    }
    resu <- optimWrapper(objectivefn=objfn,
                       initval=oldpar, gr=NULL,
                       chullformats=locchull,
                       control=list(fnscale=-1/blackbox.getOption("scalefactor"), trace=FALSE, maxit=10000)) ## returns in fittedNames space
    newpar <- fromFONKtoanyspace(resu$par, colnames(locchull$vertices))
    #if (sum(((newpar-oldpar)/localmaxrange)^2)<1e-08) break; ##  converged
    if (mean(((newpar[namesfixed]-fixedvec)/rangefixed)^2)<1e-04) break; ##  converged
    if (cost>1e10) break; ## suggests this cannot converge, for whatever reason (FR->FR: not quite clear)
    cost <- cost*10
  }
  ## resu$par is a full-dimensional vector but in locchull space
  vkrig <- tofullKrigingspace(resu$par) ## conversion to full kriging param space
  canon <- canonize(vkrig)$canonVP ## completion/conversion to canonical
  zut <- resu
  zut <- c(zut, par=canon, cost=cost)
  return(zut) ## zut$par is vector in canonical param space
} ## end profile fn prototype
