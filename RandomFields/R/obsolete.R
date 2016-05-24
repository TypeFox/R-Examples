
## Authors 
## Martin Schlather, schlather@math.uni-mannheim.de
##
##
## Copyright (C) 2015 Martin Schlather
##
## This program is free software; you can redistribute it and/or
## modify it under the terms of the GNU General Public License
## as published by the Free Software Foundation; either version 3
## of the License, or (at your option) any later version.
##
## This program is distributed in the hope that it will be useful,
## but WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
## GNU General Public License for more details.
##
## You should have received a copy of the GNU General Public License
## along with this program; if not, write to the Free Software
## Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.  



RFoldstyle <- function(old=TRUE) {
  RFoptions(general.spConform = !old,
            internal.warn_newstyle = old,
            internal.warn_oldstyle = !old)
  invisible(NULL)
}

OldDistribution <- function(distribution) {
  if (missing(distribution) || is.null(distribution) || is.na(distribution))
    return(NULL)
  names <- c("Gauss", "MaxStable")
  Distr <- list("RMgauss", NULL) ## obsolete
  nr <- pmatch(distribution, names)
  if (!is.finite(nr)) stop("unknown method")
  return(Distr[[nr]])
}

OldMethod <- function(method) {
  if (missing(method) || is.null(method) || is.na(method)) return("RPgauss")
  names <-  c("circulant embed", 
              "cutoff CE",
              "intrinsic CE",
              "TBM2", 
              "TBM3",
              "spectral TBM", 
              "direct decomp.",
              "sequential",
              "Markov",
              "average",
              "nugget", 
              "coins",
              "hyperplanes",
              ##
              "max.MPP",
              "extremalGauss",
              "BrownResnick",
              "BooleanFunction")
  Meth <- c("circulant", "cutoff", "intrinsic", "tbm", "tbm",
            "spectral", "direct",  "sequential",
            "Markov does not work anymore",
            "average", "nugget", "coins", "hyperplane",
            "smith", "schlather", "brownresnick", "schlather")
  nr <- pmatch(method, names)
  if (!is.finite(nr)) stop("unknown method")
  return(paste("RP", Meth[nr], sep=""))
}


Variogram <- function(x,
                      model, param=NULL,
              dim=if (!missing(Distances)) {if (is.matrix(x)) ncol(x) else 1},
                      ## y=NULL,
                      Distances) {
  if (RFoptions()$internal$warn_oldstyle)
    warning("The function is obsolete. Use 'RFvariogram' instead.")
  if (!missing(param) && is.na(param[1])) param[1] <- 0

##  Print(model, param)
  model <- PrepareModel(model, param)
  RFvariogram(x=x, model=model, dim=dim, distances=Distances)
}


Covariance <- CovarianceFct <-
  function(x, y=NULL,
           model, param=NULL,
           dim=if (!missing(Distances)) {if (is.matrix(x)) ncol(x) else 1},
           Distances, fctcall=c("Cov", "Variogram", "CovMatrix")) {
    if (RFoptions()$internal$warn_oldstyle)
      warning("The function is obsolete. Use 'RFcov' instead")
    if (!missing(param) && is.na(param[1])) param[1] <- 0
    model <- PrepareModel(model, param)
    fctcall <- match.arg(fctcall)
    rfeval(x=x, y=y, model=model, dim=dim, distances=Distances, fctcall=fctcall)
  }

CovMatrix <- function(x, y=NULL,
           model, param=NULL,
           dim=if (!missing(Distances)) {if (is.matrix(x)) ncol(x) else 1},
           Distances) {
  if (!missing(param) && is.na(param[1])) param[1] <- 0
  model <- PrepareModel(model, param)
  if (RFoptions()$internal$warn_oldstyle)
    warning("The function is obsolete. Use 'RFcovmatrix' instead")
  RFcovmatrix(x=x, y=y, model=model, dim=dim, distances=Distances)
}

DoSimulateRF <- function (n = 1, register = 0, paired=FALSE, trend=NULL) {
  if (!is.null(trend))
    stop("parameter trend in DoSimulateRF cannot be used anymore")
  RFoptOld <-
    internal.rfoptions(register=register, gauss.paired=paired, spConform=FALSE)
  on.exit(RFoptions(LIST=RFoptOld[[1]]))

 if (RFoptOld[[1]]$internal$warn_oldstyle)
    warning("The function is obsolete.\nSee documentation of 'RFSimulate' (advanced part) instead.")

   RFsimulate(n=n)
}


InitSimulateRF <- function (x, y = NULL, z = NULL, T=NULL,
                            grid=!missing(gridtriple), model, param,
                            trend,  method = NULL,
                            register = 0, gridtriple, distribution=NA) {
  RFoptOld <-
    internal.rfoptions(register=register, #gauss.method=method,
                       spConform=FALSE)
  on.exit(RFoptions(LIST=RFoptOld[[1]]))

  if (RFoptOld[[2]]$internal$warn_oldstyle)
    warning("This function is obsolete. Use 'RFsimulate' instead.")

  OldDistribution(distribution) ## only check whether input is oK
  
  model <- c(list(OldMethod(method)),
             list(PrepareModel(model, param, trend)))

  if (!missing(gridtriple)) {
    if (gridtriple) {
      if (!is.null(x)) {
        if (is.matrix(x)) {
          x <- apply(x, 2,
                     function(s) {c(s[1], s[3], length(seq(s[1], s[2], s[3])))})
        } else {
          x <- seq(x[1], x[2], x[3])
          if (length(y)!=0)  y <- seq(y[1], y[2], y[3])
          if (length(z)!=0)  z <- seq(z[1], z[2], z[3])
        }
      }
    }
  }
  if (length(T)!=0)  T <- c(T[1], T[3], length(seq(T[1], T[2], T[3])))

  p <- list("Simulate", PrepareModel2(model))
  rfInit(model=p, x=x, y =y, z = z, T=T, grid=grid) 
}



InitGaussRF <- function(x, y = NULL, z = NULL, T=NULL,
                        grid = !missing(gridtriple),
                        model, param,
                        trend=NULL, method = NULL,
                        register = 0, gridtriple) {
  InitSimulateRF(x=x, y=y, z=z,  T=T, 
                        grid=grid, model=model, param=param, trend=trend,
                        method=method,
                        register=register,
                        gridtriple=gridtriple, distribution="Gauss")
}


GaussRF <- function (x, y = NULL, z = NULL, T=NULL,
          grid=!missing(gridtriple), model, param, trend=NULL, method = NULL, 
          n = 1, register = 0, gridtriple,
          paired=FALSE, PrintLevel=1, Storing=TRUE, ...) { #
  if (RFoptions()$internal$warn_oldstyle)
    warning("The function is obsolete. Use 'RFsimulate' instead.")


  model <- c(list(OldMethod(method)),
             list(PrepareModel(model, param, trend)))

#  Print(model)

 
  if (!missing(gridtriple)) {
    if (gridtriple) {
      if (!is.null(x)) {
        if (is.matrix(x)) {
          x <- apply(x, 2,
                     function(s) {c(s[1], s[3], length(seq(s[1], s[2], s[3])))})
        } else {
          x <- seq(x[1], x[2], x[3])
          if (length(y)!=0)  y <- seq(y[1], y[2], y[3])
          if (length(z)!=0)  z <- seq(z[1], z[2], z[3])
        }
      }
    }
  }
  if (length(T)!=0)  T <- c(T[1], T[3], length(seq(T[1], T[2], T[3])))

  RFoptOld <-
    internal.rfoptions(register=register, gauss.paired=paired,
                      # gauss.method=method,
                       spConform=FALSE)
  on.exit(RFoptions(LIST=RFoptOld[[1]]))
 
  RFsimulate(model=model, x=x, y=y, z=z, T=T, grid=grid, n=n,
             #gauss.method=method,
             register=register, gauss.paired=paired,
             printlevel = PrintLevel, storing=Storing,
             ...)

 ##   str(RFoptions())
 
}




## it does not make sense to me at the moment that a space-time model
## for extremes is defined.

InitMaxStableRF <- function(x, y = NULL, z = NULL, grid=NULL, model, param,
                            maxstable,
                            method = NULL,
                            register = 0, gridtriple = FALSE) {
  meth <- if (!is.null(method)) method else maxstable
  if (is.null(meth)) stop("method not given")
  model <- c(list(OldMethod(meth)),
             list(PrepareModel(model, param)))

  if (!missing(gridtriple)) {
    if (gridtriple) {
      if (!is.null(x)) {
        if (is.matrix(x)) {
          x <- apply(x, 2,
                     function(s) {c(s[1], s[3], length(seq(s[1], s[2], s[3])))})
        } else {
          x <- seq(x[1], x[2], x[3])
          if (length(y)!=0)  y <- seq(y[1], y[2], y[3])
          if (length(z)!=0)  z <- seq(z[1], z[2], z[3])
         }
      }
    }
  }
  if (length(T)!=0)  T <- c(T[1], T[3], length(seq(T[1], T[2], T[3])))

  RFoptOld <-
    internal.rfoptions(register=register, #gauss.method=method,
                       spConform=FALSE)
  on.exit(RFoptions(LIST=RFoptOld[[1]]))
  
  p <- list("Simulate", PrepareModel2(model))
  rfInit(model=p, x=x, y=y, z=z, grid=grid)
}

  
MaxStableRF <- function (x, y = NULL, z = NULL, grid=NULL,
                         model, param, maxstable,
                         method = NULL,
                         n = 1, register = 0,
                         gridtriple = FALSE,
                         ...) {
  meth <- if (!is.null(method)) method else maxstable
  if (is.null(meth)) stop("method not given")
  model <- c(list(OldMethod(meth)),
             list(PrepareModel(model, param)))

  if (!missing(gridtriple)) {
    if (gridtriple) {
      if (!is.null(x)) {
        if (is.matrix(x)) {
          x <- apply(x, 2,
                     function(s) {c(s[1], s[3], length(seq(s[1], s[2], s[3])))})
        } else {
          x <- seq(x[1], x[2], x[3])
          if (length(y)!=0)  y <- seq(y[1], y[2], y[3])
          if (length(z)!=0)  z <- seq(z[1], z[2], z[3])
        }
      }
    }
  }
  
  RFoptOld <-
   if (n>1) internal.rfoptions(..., register=register, #gauss.method=method,
                               spConform=FALSE, storing=TRUE)
   else internal.rfoptions(..., register=register, #gauss.method=method,
                           spConform=FALSE)
  on.exit(RFoptions(LIST=RFoptOld[[1]]))
   
  return(RFsimulate(model=model, x=x, y=y, z=z, grid=grid, n=n))
  
}


EmpiricalVariogram <-
  function (x, y = NULL, z = NULL, T=NULL, data, grid=NULL, bin, gridtriple = FALSE,
            phi,  ## phi[1] erste richtung, phi[2] : anzahl der richtungen
            theta, ## aehnlich
            deltaT ##  deltaT[1] max abstand, deltaT[2] : gitterabstand
            )
{
  if (RFoptions()$internal$warn_oldstyle)
    warning("This function is obsolete. Use RFempiricalvariogram instead.")
  
  if (!missing(gridtriple)) {
    if (gridtriple) {
      if (!is.null(x)) {
        if (is.matrix(x)) {
          x <- apply(x, 2,
                     function(s) {c(s[1], s[3], length(seq(s[1], s[2], s[3])))})
        } else {
          x <- seq(x[1], x[2], x[3])
          if (length(y)!=0)  y <- seq(y[1], y[2], y[3])
          if (length(z)!=0)  z <- seq(z[1], z[2], z[3])
        }
      }
    }
  }
  if (length(T)!=0)  T <- c(T[1], T[3], length(seq(T[1], T[2], T[3])))
  
  RFempiricalvariogram(x=x, y=y, z=z, T=T, data=data, grid=grid,
                       bin=bin, phi=phi, theta=theta, deltaT=deltaT)
}



Kriging <- function(krige.method, x, y=NULL, z=NULL, T=NULL,
                    grid=NULL, gridtriple=FALSE,
                    model, param, given, data, trend=NULL,            
                    pch=".", return.variance=FALSE,
                    allowdistanceZero = FALSE, cholesky=FALSE) {
                    
  if (!is.null(trend)) 
    stop("in the obsolete setting, Kriging may not be used with trend!\n")                  
                    
  model <- PrepareModel(model, param, trend)
  if (!missing(gridtriple)) {
    if (gridtriple) {
      if (!is.null(x)) {
        if (is.matrix(x)) {
          x <- apply(x, 2,
                     function(s) {c(s[1], s[3], length(seq(s[1], s[2], s[3])))})
        } else {
          x <- seq(x[1], x[2], x[3])
          if (length(y)!=0)  y <- seq(y[1], y[2], y[3])
          if (length(z)!=0)  z <- seq(z[1], z[2], z[3])
        }
      }
    }
  }
  if (length(T)!=0)  T <- c(T[1], T[3], length(seq(T[1], T[2], T[3])))

  RFoptOld <- internal.rfoptions(general.pch=pch, spConform=FALSE,
                                 return_variance=return.variance,
                                 allowdistanceZero=allowdistanceZero)
  on.exit(RFoptions(LIST=RFoptOld[[1]]))
  
  data <- cbind(given, data)
  colnames(data) <- c(rep("", ncol(given)), "data")
  
  RFinterpolate(x=x, y=y, z=z, T=T, grid=grid, model=model, data=data)
}


CondSimu <- function(krige.method, x, y=NULL, z=NULL, T=NULL,
                     grid=NULL, gridtriple=FALSE,
                     model, param, method=NULL,
                     given, data, trend=NULL,
                     n=1, register=0, 
                     err.model=NULL, err.param=NULL, err.method=NULL,
                     err.register=1, 
                     tol=1E-5, pch=".", #
                     paired=FALSE,
                     na.rm=FALSE #
                     ) {
  RFoptOld <- internal.rfoptions(register=register, gauss.paired=paired,
                                        # gauss.method=method,
                                 general.errregister=err.register,
                                 spConform=FALSE)
  on.exit(RFoptions(LIST=RFoptOld[[1]]))
  
  if (RFoptOld[[2]]$internal$warn_oldstyle)
    warning("The function is obsolete.\nUse 'RFsimulate' instead.")
  if (!is.null(err.method)) warning("err.method is ignored.")
  
  ##  Print(krige.method, is.character(krige.method), (is.na(pmatch(krige.method, c("S","O")))))
  
 if (is.character(krige.method) && is.na(pmatch(krige.method, c("S","O"))))
   stop("Sorry. The parameters of the function `CondSimu' as been redefined. Use `krige.method' instead of `method'.")

  model <- PrepareModel(model, param, trend)
  if (!is.null(trend)) stop("trend cannot be given anymore.")
  err.model <- PrepareModel(err.model, err.param)
  if (err.model[[1]]=="+" && length(err.model) == 3 &&
      err.model[[3]][[1]]=="trend" && err.model[[3]]$mean != 0) {
    warning("error model has non-zero mean. Set to zero.")
    err.model <- err.model[[2]]
  } 
  if (!missing(gridtriple)) {
    if (gridtriple) {
      if (!is.null(x)) {
        if (is.matrix(x)) {
          x <- apply(x, 2,
                     function(s) {c(s[1], s[3], length(seq(s[1], s[2], s[3])))})
        } else {
          x <- seq(x[1], x[2], x[3])
          if (length(y)!=0)  y <- seq(y[1], y[2], y[3])
          if (length(z)!=0)  z <- seq(z[1], z[2], z[3])
        }
      }
    }
  }
  if (length(T)!=0)  T <- c(T[1], T[3], length(seq(T[1], T[2], T[3])))
  
 
  data <- cbind(given, data)
  if (is.null(dim(given)))
    given <- matrix(given)
  colnames(data) <- c(rep("", ncol(given)), "data")
  return(rfCondGauss(model, x=x, y=y, z=z, T=T, grid=grid, n=n,
                     data=data, err.model=err.model)$simu)
  
}


  
hurst <- function(x, y = NULL, z = NULL, data,
                  gridtriple = FALSE, sort=TRUE,
                  block.sequ = unique(round(exp(seq(log(min(3000, dim[1] / 5)),
                    log(dim[1]), len=min(100, dim[1]))))),
                  fft.m = c(1, min(1000, (fft.len - 1) / 10)),
                  fft.max.length = Inf, ## longer ts are cut down
                  method=c("dfa", "fft", "var"),
                  mode=c("plot", "interactive"),
                   pch=16, cex=0.2, cex.main=0.85,
                  PrintLevel=RFoptions()$general$printlevel,
                  height=3.5,
                  ...
                  ) {
  if (RFoptions()$internal$warn_oldstyle)
    warning("'hurst' is obsolete. Use RFhurst instead.")
  fft.len <- min(dim[1], fft.max.length)
 
  if (missing(x)) {
    gridtriple <- TRUE
    ## stopifnot(grid) 
    x <- if (is.array(data)) rbind(1, dim(data), 1) else c(1, length(data), 1)
  }
  
  if (gridtriple) {
    if (!is.null(x)) {
      if (is.matrix(x)) {
        x <- apply(x, 2,
                   function(s) {c(s[1], s[3], length(seq(s[1], s[2], s[3])))})
      } else {
        x <- seq(x[1], x[2], x[3])
        if (length(y)!=0)  y <- seq(y[1], y[2], y[3])
        if (length(z)!=0)  z <- seq(z[1], z[2], z[3])
      }
    }
  }
  
  RFhurst(x=x, y=y, z=z, data=data, sort=sort,
          block.sequ=block.sequ, fft.m=fft.m, fft.max.length=fft.max.length,
          #gauss.method=method,
          #mode=mode, pch=pch, cex=cex, cex.main=cex.main,
          printlevel=PrintLevel, height=height, ...)
}



fractal.dim <-
  function(x, y = NULL, z = NULL, data,
           grid=TRUE, gridtriple = FALSE,
           bin,
           vario.n=5,
           sort=TRUE,
           #box.sequ=unique(round(exp(seq(log(1),
           #  log(min(dim - 1, 50)), len=100)))),
           #box.enlarge.y=1,
           #range.sequ=unique(round(exp(seq(log(1),
           #  log(min(dim - 1, 50)), len=100)))),
           fft.m = c(65, 86), ## in % of range of l.lambda
           fft.max.length=Inf,
           fft.max.regr=150000,
           fft.shift = 50, # in %; 50:WOSA; 100: no overlapping
           method=c("variogram", "fft"),# "box","range", not correctly implement.
           mode=c("plot", "interactive"),
           pch=16, cex=0.2, cex.main=0.85,
           PrintLevel = RFoptions()$general$printlevel,
           height=3.5,
           ...) {

    if (RFoptions()$internal$warn_oldstyle)
      warning("`fractal.dim' is obsolete. Use RFfractaldim instead.")
    
    if (missing(x)) {
      gridtriple <- TRUE
      stopifnot(grid)
      x <- if (is.array(data)) rbind(1, dim(data), 1) else c(1, length(data), 1)
    }
    
    if (gridtriple) {
      if (!is.null(x)) {
        if (is.matrix(x)) {
          x <- apply(x, 2,
                     function(s) {c(s[1], s[3], length(seq(s[1], s[2], s[3])))})
        } else {
          x <- seq(x[1], x[2], x[3])
          if (length(y)!=0)  y <- seq(y[1], y[2], y[3])
          if (length(z)!=0)  z <- seq(z[1], z[2], z[3])
        }
      }
    }

    RFfractaldim(x=x, y=y, z=z, data=data, grid=grid, bin=bin, vario.n=vario.n,
                 sort=sort, fft.m=fft.m, fft.max.length=fft.max.length,
                 fft.max.regr=fft.max.regr, fft.shift=fft.shift,
                # gauss.method=method,
                 mode=mode, pch=pch, cex=cex, cex.main=cex.main,
                 printlevel=PrintLevel, heigth=height, ...)
    
  }




fitvario <-
  function(x, y=NULL, z=NULL, T=NULL, data, model, param,
           lower=NULL, upper=NULL, sill=NA, grid=!missing(gridtriple),
           gridtriple=FALSE,
           ...) {
    fitvario.default(x=x, y=y, z=z, T=T, data=data, model=model, param=param,
                     lower=lower, upper=upper, sill=sill,
                     grid=grid, gridtriple,     
                     ...)
  }

fitvario.default <-
  function(x, y=NULL, z=NULL, T=NULL, data, model, param,
           grid=!missing(gridtriple), gridtriple=FALSE,     
           trend = NULL,
##         BoxCox ((x+c)^l - 1) / l; log(x+c); with c==1
##         BC.lambda : NA / value
##         BC.c: nur value
           BC.lambda, ## if missing then no BoxCox-Trafo
           BC.lambdaLB=-10, BC.lambdaUB=10,
           lower=NULL, upper=NULL, sill=NA,
           use.naturalscaling=FALSE,
           ## speed=FALSE, i.e. coordniates will be saved in GATTER
           PrintLevel=RFoptions()$general$printlevel, optim.control=NULL,
           bins=20, nphi=1, ntheta=1, ntime=20,
           bin.dist.factor=0.5,
           upperbound.scale.factor=3, lowerbound.scale.factor=3, 
           lowerbound.scale.LS.factor=5,
           upperbound.var.factor=10, lowerbound.var.factor=100,     
           lowerbound.sill=1E-10, scale.max.relative.factor=1000,
           minbounddistance=0.001, minboundreldist=0.02,
           approximate.functioncalls=50, 
           minmixedvar=1/1000, maxmixedvar=1000,
           pch=RFoptions()$general$pch, 
           transform=NULL,  standard.style=NA,
     ##      var.name="X", time.name="T",
           lsq.methods=c("self", "plain", "sqrt.nr", "sd.inv", "internal"),
           ## "internal" : name should not be changed; should always be last
           ##              method!
           mle.methods=c("ml"), # "reml", "rml1"),
           cross.methods=NULL,
  #       cross.methods=c("cross.sq", "cross.abs", "cross.ign", "cross.crps"),
           users.guess=NULL,  only.users = FALSE,
           distances=NULL, truedim,
           solvesigma = NA, 
           allowdistanceZero = FALSE,
           na.rm = TRUE) { ## do not use FALSE for mixed models !

     
  RFoptOld <-
    internal.rfoptions(fit.boxcox_lb=BC.lambdaLB,
            fit.boxcox_ub=BC.lambdaUB,
     #       fit.sill=sill,
            fit.use_naturalscaling= use.naturalscaling,
            printlevel=PrintLevel, 
            fit.bins=bins,
            fit.nphi=nphi,
            fit.ntheta=ntheta,
            fit.ntime=ntime,
            fit.bin_dist_factor = bin.dist.factor,
            fit.upperbound_scale_factor = upperbound.scale.factor,
            fit.lowerbound_scale_factor = lowerbound.scale.factor,
            fit.lowerbound_scale_ls_factor=lowerbound.scale.LS.factor,
            fit.upperbound_var_factor=upperbound.var.factor,
            fit.lowerbound_var_factor=lowerbound.var.factor,     
     #       fit.lowerbound_sill=lowerbound.sill,
            fit.scale_max_relative_factor=scale.max.relative.factor,
            fit.minbounddistance=minbounddistance,
            fit.minboundreldist=minboundreldist,
            fit.approximate_functioncalls=approximate.functioncalls,
            fit.minmixedvar=minmixedvar, fit.maxmixedvar=maxmixedvar,
            pch=pch, 
    #        fit.optim_var_elimination=
    #                   if(is.na(standard.style)) 'try' else if (standard.style)
    #                   'yes' else 'never',
            ## var.name="X", time.name="T",
            fit.only_users = only.users,
    #        fit.solvesigma = solvesigma , 
            allowdistanceZero = allowdistanceZero,
                       na_rm = na.rm)
  on.exit(RFoptions(LIST=RFoptOld[[1]]))


  if (RFoptOld[[2]]$internal$warn_oldstyle)
      warning("fitvario is obsolete. Use RFfit instead.")
  model <- PrepareModel(model, param, trend)
    
  if (gridtriple) {
    if (!is.null(x)) {
      if (is.matrix(x)) {
        x <- apply(x, 2,
                   function(s) {c(s[1], s[3], length(seq(s[1], s[2], s[3])))})
      } else {
        x <- seq(x[1], x[2], x[3])
        if (length(y)!=0)  y <- seq(y[1], y[2], y[3])
        if (length(z)!=0)  z <- seq(z[1], z[2], z[3])
       }
    }
  }
  if (length(T)!=0)  T <- c(T[1], T[3], length(seq(T[1], T[2], T[3])))
 

  if (!missing(param)) {
    if (is.numeric(lower)) lower <- PrepareModel(model=model, param=lower)
    if (is.numeric(upper)) upper <- PrepareModel(model=model, param=upper)
  }
   
  RFfit(x=x, y=y, z=z, T=T, data=data, model=model,
        lower=lower, upper=upper, grid=grid,boxcox=BC.lambda,
        sub.methods=lsq.methods, methods=mle.methods,
        users.guess=users.guess, distances=distances,
        dim= truedim,
        optim.control=optim.control, transform=transform,
        spConform=FALSE)
}

DeleteAllRegisters <- function() {
  old <- RFoptions()$general$storing
  RFoptions(storing=FALSE, storing=old)
}

DeleteRegister <- function(nr=0){
 DeleteAllRegisters() 
}


RFparameters <- function(...) {
  l <- list(...)
  nl <- names(l)
  names <- c("Storing", "PrintLevel") #
  if (any(is.na(pmatch(nl, names))))
    warning("some options of 'RFparameters' in the package 'RandomFields' are not recognized anymore. Use 'RFoptions' instead")
  if ("Storing" %in% nl) RFoptions(storing=l$Storing)
  if ("PrintLevel" %in% nl) RFoptions(printlevel=l$PrintLevel)
  
}


plotWithCircles <- function(data, factor=1.0,
                            xlim=range(data[,1])+c(-maxr,maxr),
                            ylim=range(data[,2])+c(-maxr,maxr),
                            col=1, fill=0, ...) {
  ## marked point process: presents positive values of data as radii of circles
  CIRLLE.X <- cos(seq(0,2*pi,l=20))
  CIRLLE.Y <- sin(seq(0,2*pi,l=20))
  circle <- function(x,r) { polygon(x[1]+ r* CIRLLE.X,x[2]+ r* CIRLLE.Y,
                                    col=fill, border=col) }
  ##r <- x$NormedData - min(x$NormedData) +1
  ##r <- r/max(r)/nrow(x$coord) * diff(xlim) * diff(ylim) * 2.5;
  maxr <- max(data[,3])
  plot(Inf, Inf, xlim=xlim, ylim=ylim, xlab="", ylab="",...)
  for (i in 1:nrow(data)) { circle(data[i,c(1,2)], factor*data[i,3]) }
}

