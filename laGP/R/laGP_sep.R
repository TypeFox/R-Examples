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
# Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301 USA
#
# Questions? Contact Robert B. Gramacy (rbgramacy@chicagobooth.edu)
#
#*******************************************************************************



## laGPsep:
##
## C-version of sequential design loop for prediction at Xref

laGPsep <- function(Xref, start, end, X, Z, d=NULL, g=1/1000,
                 method=c("alc", "alcray", "nn"), Xi.ret=TRUE, 
                 close=min(1000*if(method == "alcray") 10 else 1, nrow(X)), 
                 numrays=ncol(X), rect=NULL, verb=0)
  {
    ## argument matching and numerifying
    method <- match.arg(method)
    if(method == "alc") imethod <- 1
    else if(method == "alcray") imethod <- 2
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
    } else { 
      if(!is.null(rect)) warning("rect only used by alcray method"); 
      rect <- 0 
    }

    ## sanity checks on input dims
    if(start < 6 || end <= start) stop("must have 6 <= start < end")
    if(ncol(Xref) != ncol(X)) stop("bad dims")
    if(length(Z) != nrow(X)) stop("bad dims")
    if(nrow(X) <= end) stop("nrow(X) <= end so nothing to do")

    ## process the d argument
    d <- darg(d, X)
    if(length(d$start) == 1) d$start <- rep(d$start, ncol(X))
    else if(length(d$start) != ncol(X)) 
      stop("d$start should be scalar or length ncol(X)")
    ## process the g argument
    g <- garg(g, Z)
    if(length(g$start) != 1) stop("g$start should be scalar")

    ## convert to doubles
    m <- ncol(X)
    dd <- c(d$start, d$mle, rep(d$min, m), rep(d$max, m), d$ab)
    dg <- c(g$start, g$mle, g$min, g$max, g$ab)

    ## sanity checks on controls
    if(!(is.logical(Xi.ret) && length(Xi.ret) == 1))
      stop("Xi.ret not a scalar logical")
    
    ## for timing
    tic <- proc.time()[3]
    
    out <- .C("laGPsep_R",
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
              numrays = as.integer(numrays),
              rect = as.double(t(rect)),
              verb = as.integer(verb),
              Xi.ret = as.integer(Xi.ret),
              Xi = integer(end*Xi.ret),
              mean = double(nref),
              s2 = double(nref),
              df = double(1),
              dmle = double(m * d$mle),
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
    if(d$mle) mle <- data.frame(d=matrix(out$dmle, nrow=1), dits=out$dits)
    if(g$mle) mle <- cbind(mle, data.frame(g=out$gmle, gits=out$gits)) 
    outp$mle <- mle
    if(Xi.ret) outp$Xi <- out$Xi + 1

    ## add ray info?
    if(method == "alcray") outp$numrays <- numrays

    ##return
    return(outp)
  }


## laGPsep.R:
##
## and R-loop version of the laGPsep function; the main reason this is 
## much slower than the C-version (laGPsep) is that it must pass/copy
## a big X-matrix each time it is called

laGPsep.R <- function(Xref, start, end, X, Z, d=NULL, g=1/1000,
                   method=c("alc", "alcray", "nn"), 
                   Xi.ret=TRUE, pall=FALSE, 
                   close=min(1000*if(method == "alcray") 10 else 1, nrow(X)),
                   parallel=c("none", "omp"), numrays=ncol(X), 
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

    ## process the d argument
    d <- darg(d, X)
    if(length(d$start) == 1) d$start <- rep(d$start, ncol(X))
    else if(length(d$start) != ncol(X)) 
      stop("d$start should be scalar or length ncol(X)")
    ## process the g argument
    g <- garg(g, Z)
    if(length(g$start) != 1) stop("g$start should be scalar")

    ## check Xi.ret argument
    if(!( is.logical(Xi.ret) && length(Xi.ret) == 1))
      stop("Xi.ret not a scalar logical")
    if(Xi.ret) Xi.ret <- rep(NA, end)
    else Xi.ret <- NULL

    ## for timing
    tic <- proc.time()[3]

    ## sorting to Xref location, and building new GPsep
    dst <- drop(distance(Xref, X))
    if(is.matrix(dst)) dst <- apply(dst, 2, min)
    cands <- order(dst)
    Xi <- cands[1:start]

    ## building a new GP with closest Xs to Xref, no derivatives
    gpsepi <- newGPsep(X[Xi,,drop=FALSE], Z[Xi], d=d$start, g=g$start) 
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
        pall[t-start,] <- predGPsep(gpsepi, Xref, lite=TRUE)

      ## calc ALC to reference
      if(method == "alcray") { 
        offset <- ((t-start) %% floor(sqrt(t-start))) + 1
        w <- lalcrayGPsep(gpsepi, Xref, X[cands,,drop=FALSE], rect, offset, numrays, verb=verb-2)
      } else {
        if(method == "alc") 
          als <- alcGPsep(gpsepi, X[cands,,drop=FALSE], Xref, 
                          parallel=parallel, verb=verb-2)
        else als <- c(1, rep(0, length(cands)-1)) ## nearest neighbor
        als[!is.finite(als)] <- NA
        w <- which.max(als)
      }

      ## add the chosen point to the GPsep fit
      updateGPsep(gpsepi, matrix(X[cands[w],], nrow=1), Z[cands[w]], verb=verb-1)
      if(!is.null(Xi.ret)) Xi.ret[t] <- cands[w]
      cands <- cands[-w]
    }

    ## maybe do post-MLE calculation 
    mle <- mleGPsep.switch(gpsepi, method, d, g, verb)

    ## Obtain final prediction
    outp <- predGPsep(gpsepi, Xref, lite=(nrow(Xref)==1))
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
    deleteGPsep(gpsepi)
    
    return(outp)
  }


## aGPsep.R:
##
## loops over all predictive locations XX and obtains adaptive approx
## kriging equations for each based on localized subsets of (X,Z); 
## the main reason this is much slower than the C-version (aGPsep) is 
## that it must pass/copy a big X-matrix each time it is called

aGPsep.R <- function(X, Z, XX, start=6, end=50, d=NULL, g=1/1000,
                  method=c("alc", "alcray", "nn"), Xi.ret=TRUE, 
                  close=min(1000*if(method == "alcray") 10 else 1, nrow(X)),
                  numrays=ncol(X), laGPsep=laGPsep.R, verb=1)
  {
    ## sanity checks
    nn <- nrow(XX)
    m <- ncol(X)
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

    ## get d and g arguments
    d <- darg(d, X)
    g <- garg(g, Z)

    ## check d$start
    ds.norep <- d$start
    if(length(d$start) == 1) 
      d$start <- matrix(rep(d$start, m), ncol=m, nrow=nn, byrow=TRUE)
    else if(length(d$start) == m) 
      d$start <- matrix(d$start, nrow=nn, byrow=TRUE)
    else if(nrow(d$start) != nn || ncol(d$start) != m)
      stop("d$start must be a scalar, or a vector of length ncol(X), or an nrow(XX) x ncol(X) matrix")
    ## check gstart
    if(length(g$start) > 1 && length(g$start) != nn) 
      stop("g$start must be a scalar or a vector of length nrow(XX)")
    gs.norep <- g$start
    if(length(g$start) != nrow(XX)) g$start <- rep(g$start, nrow(XX))

    ## check mle
    if(d$mle) {
      dits <- ZZ.var 
      dmle <- matrix(NA, nrow=nrow(XX), ncol=ncol(X))
    } else dits <- dmle <- NULL
    if(g$mle) gits <- gmle <- ZZ.var
    else gits <- gmle <- NULL

    ## for timing
    tic <- proc.time()[3]
    
    ## now do copies and local updates for each reference location
    for(i in 1:N) {

      ## local calculation, (add/remove .R in laGP.R for R/C version)
      di <- list(start=d$start[i,], mle=d$mle, min=d$min, max=d$max, ab=d$ab)
      gi <- list(start=g$start[i], mle=g$mle, min=g$min, max=g$max, ab=g$ab)
      outp <- laGPsep(XX[i,,drop=FALSE], start, end, X, Z, d=di, g=gi, 
        method=method, Xi.ret=Xi.ret, close=close, numrays=numrays, 
        rect=rect, verb=verb-1)

      ## save MLE outputs and update gpi to use new dmle
      if(!is.null(dmle)) { dmle[i,] <- as.numeric(outp$mle[1:ncol(X)]); dits[i] <- outp$mle$dits }
      if(!is.null(gmle)) { gmle[i] <- outp$mle$g; gits[i] <- outp$mle$gits }

      ## extract predictive equations
      ZZ.mean[i] <- outp$mean
      ZZ.var[i] <- outp$s2 * outp$df / (outp$df-2)

      ## save Xi; Xi.ret checked in laGP.R
      if(Xi.ret) Xi[i,] <- outp$Xi

      ## print progress
      if(verb > 0) {
        cat("i = ", i, " (of ", N, ")", sep="")
        if(d$mle) cat(", d = (", paste(signif(dmle[i,], 5), collapse=", "), "), its = ", dits[i], sep="")
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


## aGPsep:
##
## using C: loops over all predictive locations XX and obtains adaptive
## approx kriging equations for each based on localized subsets of (X,Z)

aGPsep <- function(X, Z, XX, start=6, end=50, d=NULL, g=1/1000,
                method=c("alc", "alcray", "nn"), Xi.ret=TRUE, 
                close=min(1000*if(method == "alcray") 10 else 1, nrow(X)), 
                numrays=ncol(X), omp.threads=1, verb=1)
  {
    ## sanity checks
    nn <- nrow(XX)
    m <- ncol(X)
    if(ncol(XX) != m) stop("mismatch XX and X cols")
    if(nrow(X) != length(Z)) stop("length(Z) != nrow(X)")
    if(end-start <= 0) stop("nothing to do")

    ## numerify method
    method <- match.arg(method)
    if(method == "alc") imethod <- 1
    else if(method == "alcray") imethod <- 2
    else imethod <- 5

    ## calculate rectangle if using alcray
    if(method == "alcray") {
      rect <- apply(X, 2, range)
      if(nrow(rect) != 2 || ncol(rect) != m)
        stop("bad rect dimensions, must be 2 x ncol(X)")
      if(length(numrays) != 1 || numrays < 1)
        stop("numrays should be an integer scalar >= 1")
    } else rect <- 0

    ## check Xi.ret argument
    if(!(is.logical(Xi.ret) && length(Xi.ret) == 1))
      stop("Xi.ret not a scalar logical")

        ## get d and g arguments
    d <- darg(d, X)
    dd <- c(d$mle, rep(d$min, m), rep(d$max, m), d$ab)
    g <- garg(g, Z)
    dg <- c(g$mle, g$min, g$max, g$ab)

    ## check d$start
    ds.norep <- d$start
    if(length(d$start) == 1) 
      d$start <- matrix(rep(d$start, m), ncol=m, nrow=nn, byrow=TRUE)
    else if(length(d$start) == m) 
      d$start <- matrix(d$start, nrow=nn, byrow=TRUE)
    else if(nrow(d$start) != nn || ncol(d$start) != m)
      stop("d$start must be a scalar, or a vector of length ncol(X), or an nrow(XX) x ncol(X) matrix")
    ## check gstart
    if(length(g$start) > 1 && length(g$start) != nn) 
      stop("g$start must be a scalar or a vector of length nrow(XX)")
    gs.norep <- g$start
    if(length(g$start) != nrow(XX)) g$start <- rep(g$start, nrow(XX))

    ## check OMP argument
    if(length(omp.threads) != 1 || omp.threads < 1)
      stop("omp.threads should be a positive scalar integer")

    ## for timing
    tic <- proc.time()[3]

    ## calculate the kriging equations separately
    out <- .C("aGPsep_R",
              m = as.integer(m),
              start = as.integer(start),
              end = as.integer(end),
              XX = as.double(t(XX)),
              nn = as.integer(nn),
              n = as.integer(nrow(X)),
              X = as.double(t(X)),
              Z = as.double(Z),
              dstart = as.double(t(d$start)),
              darg = as.double(dd),
              g = as.double(g$start),
              garg = as.double(dg),
              imethod = as.integer(imethod),
              close = as.integer(close),
              omp.threads = as.integer(omp.threads),
              numrays = as.integer(numrays),
              rect = as.double(t(rect)),
              verb = as.integer(verb),
              Xi.ret = as.integer(Xi.ret),
              Xi = integer(end*Xi.ret*nn),
              mean = double(nn),
              var = double(nn),
              dmle = double(nn * d$mle * m),
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
    if(d$mle) {
      outp$mle <- data.frame(d=matrix(out$dmle, ncol=m, byrow=TRUE), 
                             dits=out$dits)
    }
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
