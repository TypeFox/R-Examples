#*******************************************************************************
#
# Local Approximate Gaussian Process Regression
# Copyright (C) 2013, The University of Chicago
#
# This library is free software; you can redistribute it and/or
# modify it under the terms of the GNU Lesser General Public
# License as published by the Free Software Foundation; either
# version 2.1 of the License, or (at your option) any later version.
#
# This library is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
# Lesser General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public
# License along with this library; if not, write to the Free Software
# Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
#
# Questions? Contact Robert B. Gramacy (rbgramacy@chicagobooth.edu)
#
#*******************************************************************************


## laGP:
##
## C-version of sequential design loop for prediction at Xref

laGP <- function(Xref, start, end, X, Z, d=NULL, g=1/1000,
                 method=c("alc", "alcray", "mspe", "nn", "efi"), Xi.ret=TRUE, 
                 close=min(1000*if(method == "alcray") 10 else 1, nrow(X)), 
                 alc.gpu=FALSE, numrays=ncol(X), rect=NULL, verb=0)
  {
    ## argument matching and numerifying
    method <- match.arg(method)
    if(method == "alc") imethod <- 1
    else if(method == "alcray") imethod <- 2
    else if(method == "mspe") imethod <- 3
    else if(method == "efi") imethod <- 4
    else imethod <- 5

    ## massage Xref
    if(!is.matrix(Xref)) Xref <- matrix(Xref, ncol=ncol(X))
    nref <- nrow(Xref)

    ## calculate rectangle if using alcray
    if(method == "alcray") {
      if(is.null(rect)) rect <- apply(X, 2, range)
      if(nrow(Xref) != 1) stop("alcray only implemented for nrow(Xref) = 1")
      if(nrow(rect) != 2 || ncol(rect) != ncol(X))
        stop("bad rect dimensions, must be 2 x ncol(X)")
      if(length(numrays) != 1 || numrays < 1)
        stop("numrays should be an integer scalar >= 1")
    } else { if(!is.null(rect)) warning("rect only used by alcray method"); rect <- 0 }

    ## sanity checks on input dims
    if(start < 6 || end <= start) stop("must have 6 <= start < end")
    if(ncol(Xref) != ncol(X)) stop("bad dims")
    if(length(Z) != nrow(X)) stop("bad dims")
    if(nrow(X) <= end) stop("nrow(X) <= end so nothing to do")

    ## process the d argument
    d <- darg(d, X)
    dd <- c(d$start, d$mle, d$min, d$max, d$ab)
    g <- garg(g, Z)
    dg <- c(g$start, g$mle, g$min, g$max, g$ab)

    ## sanity checks on controls
    if(!(is.logical(Xi.ret) && length(Xi.ret) == 1))
      stop("Xi.ret not a scalar logical")
    if(length(alc.gpu) > 1 || alc.gpu < 0)
      stop("alc.gpu should be a scalar logical or scalar non-negative integer")
    
    ## for timing
    tic <- proc.time()[3]
    
    out <- .C("laGP_R",
              m = as.integer(ncol(Xref)),
              start = as.integer(start),
              end = as.integer(end),
              Xref = as.double(t(Xref)),
              nref = as.integer(nref),
              n = as.integer(nrow(X)),
              X = as.double(t(X)),
              Z = as.double(Z),
              d = as.double(dd),
              g = as.double(dg),
              imethod = as.integer(imethod),
              close = as.integer(close),
              alc.gpu = as.integer(alc.gpu),
              numrays = as.integer(numrays),
              rect = as.double(t(rect)),
              verb = as.integer(verb),
              Xi.ret = as.integer(Xi.ret),
              Xi = integer(end*Xi.ret),
              mean = double(nref),
              s2 = double(nref),
              df = double(1),
              dmle = double(1 * d$mle),
              dits = integer(1 * d$mle),
              gmle = double(1 * g$mle),
              gits = integer(1 * g$mle),
              llik = double(1),
              PACKAGE = "laGP")

    ## put timing in
    toc <- proc.time()[3]

    ## assemble output and return
    outp <- list(mean=out$mean, s2=out$s2, df=out$df, llik=out$llik,
                 time=toc-tic, method=method, d=d, g=g, close=close)

    ## possibly add mle and Xi info
    mle <- NULL
    if(d$mle) mle <- data.frame(d=out$dmle, dits=out$dits)
    if(g$mle) mle <- cbind(mle, data.frame(g=out$gmle, gits=out$gits)) 
    outp$mle <- mle
    if(Xi.ret) outp$Xi <- out$Xi + 1

    ## add ray info?
    if(method == "alcray") outp$numrays <- numrays

    ##return
    return(outp)
  }


## laGP.R:
##
## and R-loop version of the laGP function; the main reason this is 
## much slower than the C-version (laGP) is that it must pass/copy
## a big X-matrix each time it is called

laGP.R <- function(Xref, start, end, X, Z, d=NULL, g=1/1000,
                   method=c("alc", "alcray", "mspe", "nn", "efi"), 
                   Xi.ret=TRUE, pall=FALSE, 
                   close=min(1000*if(method == "alcray") 10 else 1, nrow(X)),
                   parallel=c("none", "omp", "gpu"), numrays=ncol(X), 
                   rect=NULL, verb=0)
  {
    ## argument matching
    method <- match.arg(method)
    parallel <- match.arg(parallel)

    ## massage Xref
    if(!is.matrix(Xref)) Xref <- matrix(Xref, ncol=ncol(X))
    
    ## sanity checks
    if(start < 6 || end <= start) stop("must have 6 <= start < end")
    if(ncol(Xref) != ncol(X)) stop("bad dims")
    if(length(Z) != nrow(X)) stop("bad dims")
    if(nrow(X) <= end) stop("nrow(X) <= end so nothing to do")

    ## calculate rectangle if using alcray
    if(method == "alcray") {
      if(is.null(rect)) rect <- apply(X, 2, range)
      if(nrow(Xref) != 1) stop("alcray only implemented for nrow(Xref) = 1")
      if(nrow(rect) != 2 || ncol(rect) != ncol(X))
        stop("bad rect dimensions, must be 2 x ncol(X)")
      if(length(numrays) != 1 || numrays < 1)
        stop("numrays should be an integer scalar >= 1")
    } else if(!is.null(rect)) warning("rect only used by alcray method")

    ## process the d and g arguments
    d <- darg(d, X)
    g <- garg(g, Z)
    ## sanity check
    if(length(d$start) != 1 || length(g$start) != 1)
      stop("laGP starting values should be scalars")

    ## check Xi.ret argument
    if(!( is.logical(Xi.ret) && length(Xi.ret) == 1))
      stop("Xi.ret not a scalar logical")
    if(Xi.ret) Xi.ret <- rep(NA, end)
    else Xi.ret <- NULL

    ## for timing
    tic <- proc.time()[3]

    ## sorting to Xref location
    dst <- drop(distance(Xref, X))
    if(is.matrix(dst)) dst <- apply(dst, 2, min)
    cands <- order(dst)
    Xi <- cands[1:start]

    ## building a new GP with closest Xs to Xref
    gpi <- newGP(X[Xi,,drop=FALSE], Z[Xi], d=d$start, g=g$start, 
                 dK=(method != "alc" && method != "alcray"))

    ## for the output object
    if(!is.null(Xi.ret)) Xi.ret[1:start] <- Xi

    ## if pall, then predict after every iteration
    ## ONLY AVAILABLE IN THE R VERSION
    if(pall) {
      nav <- rep(NA, end-start)
      pall <- data.frame(mean=nav, s2=nav, df=nav, llik=nav)
    } else pall <- NULL

    ## determine remaining candidates
    if(close >= nrow(X)) close <- 0
    if(close > 0) {
      if(close >= nrow(X)-start)
        stop("close not less than remaining cands")
      cands <- cands[(start+1):close]
    } else cands <- cands[-(1:start)]
    
    ## set up the start and end times
    for(t in (start+1):end) {

      ## if pall then predict after each iteration
      if(!is.null(pall)) 
        pall[t-start,] <- predGP(gpi, Xref, lite=TRUE)

      ## calc ALC to reference
      if(method == "alcray") {
        offset <- ((t-start) %% floor(sqrt(t-start))) + 1
        w <- lalcrayGP(gpi, Xref, X[cands,,drop=FALSE], rect, offset, numrays, 
              verb=verb-2)
      } else {
        if(method == "alc") 
          als <- alcGP(gpi, X[cands,,drop=FALSE], Xref, parallel=parallel, 
                       verb=verb-2)
        else if(method == "mspe") 
          als <- 0.0 - mspeGP(gpi, X[cands,,drop=FALSE], Xref, verb=verb-2)
        else if(method == "efi") als <- efiGP(gpi, X[cands,,drop=FALSE])
        else als <- c(1, rep(0, length(cands)-1)) ## nearest neighbor
        als[!is.finite(als)] <- NA
        w <- which.max(als)
      }

      ## add the chosen point to the GP fit
      updateGP(gpi, matrix(X[cands[w],], nrow=1), Z[cands[w]], verb=verb-1)
      if(!is.null(Xi.ret)) Xi.ret[t] <- cands[w]
      cands <- cands[-w]
    }

    ## maybe do post-MLE calculation 
    mle <- mleGP.switch(gpi, method, d, g, verb)

    ## Obtain final prediction
    outp <- predGP(gpi, Xref, lite=(nrow(Xref)==1))
    if(!is.null(pall)) outp <- as.list(rbind(pall, outp))

    ## put timing and X info in
    toc <- proc.time()[3]
    outp$time <- toc - tic
    outp$Xi <- Xi.ret
    outp$method <- method
    outp$close <- close

    ## assign d & g
    outp$d <- d
    ## assign g
    outp$g <- g
    ## assign mle
    outp$mle <- mle

    ## add ray info?
    if(method == "alcray") outp$numrays <- numrays

    ## clean up
    deleteGP(gpi)
    
    return(outp)
  }


## aGP.R:
##
## loops over all predictive locations XX and obtains adaptive approx
## kriging equations for each based on localized subsets of (X,Z); 
## the main reason this is much slower than the C-version (aGPsep) is 
## that it must pass/copy a big X-matrix each time it is called

aGP.R <- function(X, Z, XX, start=6, end=50, d=NULL, g=1/1000,
                  method=c("alc", "alcray", "mspe", "nn", "efi"), Xi.ret=TRUE, 
                  close=min(1000*if(method == "alcray") 10 else 1, nrow(X)),
                  numrays=ncol(X), laGP=laGP.R, verb=1)
  {
    ## sanity checks
    nn <- nrow(XX)
    if(ncol(XX) != ncol(X)) stop("mismatch XX and X cols")
    if(nrow(X) != length(Z)) stop("length(Z) != nrow(X)")
    if(end-start <= 0) stop("nothing to do")

    ## check method argument
    method <- match.arg(method)

    ## calculate rectangle if using alcray
    if(method == "alcray") {
      rect <- apply(X, 2, range)
      if(nrow(rect) != 2 || ncol(rect) != ncol(X))
        stop("bad rect dimensions, must be 2 x ncol(X)")
      if(length(numrays) != 1 || numrays < 1)
          stop("numrays should be an integer scalar >= 1")
    } else rect <- NULL

    ## memory for each set of approx kriging equations
    ZZ.var <- ZZ.mean <- rep(NA, nrow(XX))

    ## other args checked in laGP.R; allocate Xi space (?)
    N <- length(ZZ.mean)
    if(Xi.ret) Xi <- matrix(NA, nrow=N, ncol=end)
    else Xi <- NULL

    ## get d argument
    d <- darg(d, X)
    g <- garg(g, Z)

    ## check d$start
    if(length(d$start) > 1 && length(d$start) != nrow(XX)) 
      stop("d$start must be a scalar or a vector of length nrow(XX)")
    ds.norep <- d$start
    if(length(d$start) != nrow(XX)) d$start <- rep(d$start, nrow(XX))
    ## check gstart
    if(length(g$start) > 1 && length(g$start) != nrow(XX)) 
      stop("g$start must be a scalar or a vector of length nrow(XX)")
    gs.norep <- g$start
    if(length(g$start) != nrow(XX)) g$start <- rep(g$start, nrow(XX))

    ## check mle
    if(d$mle) dits <- dmle <- ZZ.var
    else dits <- dmle <- NULL
    if(g$mle) gits <- gmle <- ZZ.var
    else gits <- gmle <- NULL

    ## for timing
    tic <- proc.time()[3]
    
    ## now do copies and local updates for each reference location
    for(i in 1:N) {

      ## local calculation, (add/remove .R in laGP.R for R/C version)
      di <- list(start=d$start[i], mle=d$mle, min=d$min, max=d$max, ab=d$ab)
      gi <- list(start=g$start[i], mle=g$mle, min=g$min, max=g$max, ab=g$ab)
      outp <- laGP(XX[i,,drop=FALSE], start, end, X, Z, d=di, g=gi, 
                method=method, Xi.ret=Xi.ret, close=close, numrays=numrays,
                rect=rect, verb=verb-1)

      ## save MLE outputs and update gpi to use new dmle
      if(!is.null(dmle)) { dmle[i] <- outp$mle$d; dits[i] <- outp$mle$dits }
      if(!is.null(gmle)) { gmle[i] <- outp$mle$g; gits[i] <- outp$mle$gits }

      ## extract predictive equations
      ZZ.mean[i] <- outp$mean
      ZZ.var[i] <- outp$s2 * outp$df / (outp$df-2)

      ## save Xi; Xi.ret checked in laGP.R
      if(Xi.ret) Xi[i,] <- outp$Xi

      ## print progress
      if(verb > 0) {
        cat("i = ", i, " (of ", N, ")", sep="")
        if(d$mle) cat(", d = ", dmle[i], ", its = ", dits[i], sep="")
        if(g$mle) cat(", g = ", gmle[i], ", its = ", gits[i], sep="")
        cat("\n", sep="")
      }
    }

    ## for timing
    toc <- proc.time()[3]

    ## assemble output
    d$start <- ds.norep
    g$start <- gs.norep
    r <- list(Xi=Xi, mean=ZZ.mean, var=ZZ.var, d=d, g=g, 
              time=toc-tic, method=method, close=close)
    ## add mle info?
    mle <- NULL
    if(d$mle) mle <- data.frame(d=dmle, dits=dits)
    if(g$mle) mle <- cbind(mle, data.frame(g=gmle, gits=gits))
    r$mle <- mle
    ## add ray info?
    if(method == "alcray") r$numrays <- numrays

    ## done
    return(r)
  }


## aGP:
##
## using C: loops over all predictive locations XX and obtains adaptive
## approx kriging equations for each based on localized subsets of (X,Z)

aGP <- function(X, Z, XX, start=6, end=50, d=NULL, g=1/1000,
                method=c("alc", "alcray", "mspe", "nn", "efi"), Xi.ret=TRUE, 
                close=min(1000*if(method == "alcray") 10 else 1, nrow(X)), 
                numrays=ncol(X), num.gpus=0, gpu.threads=num.gpus,
                omp.threads=if(num.gpus > 0) 0 else 1, 
		            nn.gpu=if(num.gpus > 0) nrow(XX) else 0, verb=1)
  {
    ## sanity checks
    nn <- nrow(XX)
    if(ncol(XX) != ncol(X)) stop("mismatch XX and X cols")
    if(nrow(X) != length(Z)) stop("length(Z) != nrow(X)")
    if(end-start <= 0) stop("nothing to do")

    ## numerify method
    method <- match.arg(method)
    if(method == "alc") imethod <- 1
    else if(method == "alcray") imethod <- 2
    else if(method == "mspe") imethod <- 3
    else if(method == "efi") imethod <- 4
    else imethod <- 5

    ## calculate rectangle if using alcray
    if(method == "alcray") {
      rect <- apply(X, 2, range)
      if(nrow(rect) != 2 || ncol(rect) != ncol(X))
        stop("bad rect dimensions, must be 2 x ncol(X)")
      if(length(numrays) != 1 || numrays < 1)
        stop("numrays should be an integer scalar >= 1")
    } else rect <- 0

    ## check Xi.ret argument
    if(!(is.logical(Xi.ret) && length(Xi.ret) == 1))
      stop("Xi.ret not a scalar logical")

    ## process d argument
    d <- darg(d, X)
    dd <- c(d$mle, d$min, d$max, d$ab)
    g <- garg(g, Z)
    dg <- c(g$mle, g$min, g$max, g$ab)

    ## check d$start
    if(length(d$start) > 1 && length(d$start) != nrow(XX)) 
      stop("d$start must be a scalar or a vector of length nrow(XX)")
    ds.norep <- d$start
    if(length(d$start) != nrow(XX)) d$start <- rep(d$start, nrow(XX))

    ## check g$start
    if(length(g$start) > 1 && length(g$start) != nrow(XX)) 
      stop("d$start must be a scalar or a vector of length nrow(XX)")
    gs.norep <- g$start
    if(length(g$start) != nrow(XX)) g$start <- rep(g$start, nrow(XX))

    ## check OMP argument
    if(length(omp.threads) != 1 || omp.threads < 0)
      stop("omp.threads should be a non-negative scalar integer")

    ## check gpu argument
    num.gpus <- as.integer(num.gpus)    
    if(length(num.gpus) > 1 || num.gpus < 0)
      stop("num.gpus should be a non-negative scalar integer")
    gpu.threads <- as.integer(gpu.threads)    
    if(length(gpu.threads) > 1 || gpu.threads < 0)
      stop("gpu.threads should be a non-negative scalar integer")
    if(gpu.threads < num.gpus)
      cat("NOTICE: gpu.threads < num.gpus, setting gpu.threads=num.gpus\n")
    if(num.gpus > 0 && (gpu.threads %% num.gpus != 0))
      warning("suggest gpu.threads be a multiple of num.gpus")

    ## check total threads
    if(omp.threads + gpu.threads <= 0)
      stop("must specify a positive value for one of omp.threads or gpu.threads")

    ## check nn.gpu
    if(length(nn.gpu) != 1 || nn.gpu < 0) stop("nn.gpu must be a non-negative scalar integer")
    nn.gpu <- as.integer(nn.gpu)
    if(gpu.threads == 0 && nn.gpu > 0) stop("must have nn.gpu = 0 if no GPU threads")
    if(nn.gpu > 0 && nn.gpu < nn && omp.threads == 0) 
      stop("if nn.gpu != nrow(XX) then must have omp.threads > 0")
    if(nn.gpu == nn && omp.threads > 0) {
      warning("specify nn.gpu < nrow(XX) when omp.threads > 0; setting omp.threads=0")
      omp.threads <- 0
    }

    ## for timing
    tic <- proc.time()[3]

    ## calculate the kriging equations separately
    out <- .C("aGP_R",
              m = as.integer(ncol(XX)),
              start = as.integer(start),
              end = as.integer(end),
              XX = as.double(t(XX)),
              nn = as.integer(nn),
              n = as.integer(nrow(X)),
              X = as.double(t(X)),
              Z = as.double(Z),
              d = as.double(d$start),
              darg = as.double(dd),
              g = as.double(g$start),
              garg = as.double(dg),
              imethod = as.integer(imethod),
              close = as.integer(close),
              omp.threads = as.integer(omp.threads),
              num.gpus = as.integer(num.gpus),
              gpu.threads = as.integer(gpu.threads),
              nn.gpu = as.integer(nn.gpu),
              numrays = as.integer(numrays),
              rect = as.double(t(rect)),
              verb = as.integer(verb),
              Xi.ret = as.integer(Xi.ret),
              Xi = integer(end*Xi.ret*nn),
              mean = double(nn),
              var = double(nn),
              dmle = double(nn * d$mle),
              dits = integer(nn * d$mle),
              gmle = double(nn * g$mle),
              gits = integer(nn * g$mle),
              llik = double(nn),
              PACKAGE = "laGP")
    
    ## for timing
    toc <- proc.time()[3]
 
    ## all done, return
    d$start <- ds.norep
    g$start <- gs.norep
    outp <- list(mean=out$mean, var=out$var, llik=out$llik, d=d, g=g,
                 time=toc-tic, method=method, close=close)

    ## copy MLE outputs
    outp$mle <- NULL
    if(d$mle) outp$mle <- data.frame(d=out$dmle, dits=out$dits)
    if(g$mle) {
      if(d$mle) outp$mle <- cbind(outp$mle, data.frame(g=out$gmle, gits=out$gits))
      else outp$mle <- data.frame(g=out$gmle, gits=out$gits)
    }

    ## add ray info?
    if(method == "alcray") outp$numrays <- numrays

    ## copy XI
    if(Xi.ret) outp$Xi <- matrix(out$Xi+1, nrow=nn, byrow=TRUE)
    
    return(outp)
}
