## TO DO: roll back into bbmle?
## allow multiple 'transects'?
## (i.e. if two sets of parameters are given ...)
## * robustification
##  print method
## allow manual override of ranges
## allow log-scaling?
## preserve order of parameters in 1D plots

## substitute values of parameters into full parameter vector
mkpar <- function(params,p,i) {
  params[i] <- p
  params
}

## get reasonable range for slice
## document what is done here
## implement upper bound
## more robust approaches;
## try not to assume anything about signs of parameters
## inherit bounds from fitted value
get_trange <- function(pars,  ## baseline parameter values
                       i,     ## focal parameter
                       fun,   ## objective function
                       lower=-Inf, ## lower bound 
                       upper=Inf,  ## upper bound
                       cutoff=10,  ## increase above min z-value
                       maxit=200,  ## max number of iterations
                       steptype=c("mult","addprop"),
                       step=0.1) {
    ## step possibilities: multiplicative
    ## additive (absolute scale)     [not yet implemented]
    addabs <- NULL ## fix false positive test
    steptype <- match.arg(steptype)
    v <- v0 <- fun(pars)
    lowval <- pars[i]
    it <- 1
    if (steptype=="addprop") step <- step*pars[i]
    while (it<maxit && lowval>lower && v<(v0+cutoff)) {
        lowval <- switch(steptype,
                         addabs,
                         addpropn=lowval-step,
                         mult=lowval*(1-step))
        v <- fun(mkpar(pars,lowval,i))
        it <- it+1
    }
    lowdev <- v
    lowit <- it
    upval <- pars[i]
    it <- 1
    v <- v0 <- fun(pars)
    if (upval==0) upval <- 1e-4
    while (it<maxit && v<(v0+cutoff)) {
        upval <- switch(steptype,
                        addabs,
                        addpropn=lowval+step,
                         mult=lowval*(1+step))
        v <- fun(mkpar(pars,upval,i))
        ## cat(it,upper,v,"\n")
        it <- it+1
    }
    updev <- v
    upit <- it
    c(low_val=lowval,up_valr=upval,low_dev=lowdev,up_dev=updev,
      low_it=lowit,up_it=upit)
  }

get_all_trange <- function(params,fun,lower,upper,cutoff=10,...) {
  arglist <- c(list(pars=params,fun=fun,cutoff=cutoff),list(...))
  tranges <- t(mapply(FUN=get_trange,
                      seq(length(params)),
                      lower,
                      upper,
                      MoreArgs=arglist,
                      SIMPLIFY=TRUE))
  rownames(tranges) <- names(params)
  tranges
}

## generic function (S3)
slice <- function (x, dim=1, ...)  {
    UseMethod("slice")
}

slice.mle2 <- function(x,  ...) {
    ff <- x@minuslogl
    ## vectorize objective function: return minus log L
    ## (switching to deviance screws things up???)
    ff2 <- function(p) {
        do.call(ff,as.list(p))
    }
    slice0(coef(x),ff2, ...)
}

slice0 <- function(params,fun,dim=1,params2=NULL,...) {
  if (dim==1) {
    if (is.null(params2)) {
      slice1D(params,fun,...)
    } else {
      slicetrans(params,params2,fun,...)
    }
  } else {
    if (!is.null(params2)) stop("can't do transect in 2D")
    slice2D(params,fun,...)
  }
}


plot.slice <- function(x,...) {
  switch(x$dim,xyplot(x,...),splom(x,...))
}

slicetrans <- function(params, params2, fun, extend=0.1, nt=401,
                       lower=-Inf, upper=Inf) {
  ## make sure 0/1 are included in parameter vector ...
  np <- length(params)
  extend <- rep(extend,length.out=2)
  lower <- rep(lower,length.out=np)
  upper <- rep(upper,length.out=np)
  slicep <- sort(unique(c(0,1,seq(-extend[1],1+extend[2], length=nt))))
  slicepars <- t(sapply(slicep, function(x) (1 - x) * params + x * params2))
  OK <- apply(slicepars,1,function(x) all(x>=lower & x<=upper))
  if (any(!OK)) {
    warning("some parameter sets outside of bounds were removed")
    slicep <- slicep[OK]
    slicepars <- slicepars[OK,]
  }
  v <- apply(slicepars, 1, fun)
  slices <- list(data.frame(var1="trans",x=slicep,z=v))
  r <- list(slices=slices,params=params,params2=params2,dim=1)
  class(r) <- "slice"
  r
}

  
slice1D <- function(params,fun,nt=101,
                    lower=-Inf,
                    upper=Inf,
                    verbose=TRUE,
                    tranges=NULL,
                    ...) {
  npv <- length(params)
  if (is.null(pn <- names(params))) pn <- seq(npv)
  if (is.null(tranges)) {
      tranges <- get_all_trange(params,fun,
                                rep(lower,length.out=npv),
                                rep(upper,length.out=npv),
                                ...)
  }
  slices <- vector("list",npv)
  for (i in 1:npv) {
      tvec <- seq(tranges[i,1],tranges[i,2],length=nt)
      if (verbose) cat(pn[i],"\n")
      vtmp <- sapply(tvec,
                     function(t) {
                         fun(mkpar(params,t,i))})
      slices[[i]] <- data.frame(var1=pn[i],x=tvec,z=vtmp)
  }
  r <- list(slices=slices,ranges=tranges,params=params,dim=1)
  class(r) <- "slice"
  r
}

## OLD slice method
## should probably roll this in as an option to profile
## include attribute, warning? draw differently (leave off
## conf. limit lines)
## slice <- function(fitted, ...) UseMethod("slice")

## 1D slicing implemented as in profile
sliceOld <- function (fitted, which = 1:p, maxsteps = 100,
                       alpha = 0.01, zmax = sqrt(qchisq(1 - alpha/2, p)),
                       del = zmax/5, trace = FALSE,
                       tol.newmin=0.001, ...)
{
    onestep <- function(step)
    {
        bi <- B0[i] + sgn * step * del * std.err[i]
        fix <- list(bi)
        names(fix) <- p.i
        call$fixed <- c(fix,eval(call$fixed))
        call$eval.only = TRUE
        pfit <- try(eval(call), silent=TRUE) ##
        if(inherits(pfit, "try-error")) return(NA)
        else {
            zz <- 2*(pfit@min - fitted@min)
            ri <- pv0
            ri[, names(pfit@coef)] <- pfit@coef
            ri[, p.i] <- bi
            if (zz > -tol.newmin)
                zz <- max(zz, 0)
            else stop("profiling has found a better solution, so original fit had not converged")
            z <- sgn * sqrt(zz)
            pvi <<- rbind(pvi, ri)
            zi <<- c(zi, z) ## NB global set!
        }
        if (trace) cat(bi, z, "\n")
        z
      }
    ## Profile the likelihood around its maximum
    ## Based on profile.glm in MASS
    summ <- summary(fitted)
    std.err <- summ@coef[, "Std. Error"]
    Pnames <- names(B0 <- fitted@coef)
    pv0 <- t(as.matrix(B0))
    p <- length(Pnames)
    prof <- vector("list", length = length(which))
    names(prof) <- Pnames[which]
    call <- fitted@call
    call$minuslogl <- fitted@minuslogl
    for (i in which) {
        zi <- 0
        pvi <- pv0
        p.i <- Pnames[i]
        for (sgn in c(-1, 1)) {
          if (trace)
            cat("\nParameter:", p.i, c("down", "up")[(sgn + 1)/2 + 1], "\n")
          step <- 0
          z <- 0
          ## This logic was a bit frail in some cases with
          ## high parameter curvature. We should probably at least
          ## do something about cases where the mle2 call fails
          ## because the parameter gets stepped outside the domain.
          ## (We now have.)
          call$start <- as.list(B0)
          lastz <- 0
          while ((step <- step + 1) < maxsteps && abs(z) < zmax) {
            z <- onestep(step)
            if(is.na(z)) break
            lastz <- z
          }
          if(abs(lastz) < zmax) {
            ## now let's try a bit harder if we came up short
            for(dstep in c(0.2, 0.4, 0.6, 0.8, 0.9)) {
              z <- onestep(step - 1 + dstep)
              if(is.na(z) || abs(z) > zmax) break
            }
          } else if(length(zi) < 5) { # try smaller steps
            mxstep <- step - 1
            step <- 0.5
            while ((step <- step + 1) < mxstep) onestep(step)
          }
        }
        si <- order(pvi[, i])
        prof[[p.i]] <- data.frame(z = zi[si])
        prof[[p.i]]$par.vals <- pvi[si,, drop=FALSE]
    }
    list(profile = prof, summary = summ)
}



## * is it possible to set up the 2D vectors so they include
##   the baseline value? maybe not easily ...
slice2D <- function(params,
                    fun,
                    nt=31,
                    lower=-Inf,
                    upper=Inf,
                    cutoff=10,
                    verbose=TRUE,
                    tranges=NULL,
                    ...) {
  npv <- length(params)
  if (is.null(pn <- names(params))) pn <- seq(npv)
  if (is.null(tranges)) {
      tranges <- get_all_trange(params,fun,
                                rep(lower,length.out=npv),
                                rep(upper,length.out=npv),
                                cutoff=cutoff,
                                ...)
  }
  slices <- list()
  for (i in 1:(npv-1)) {
    slices[[i]] <- vector("list",npv)
    for (j in (i+1):npv) {
      if (verbose) cat("param",i,j,"\n")
      t1vec <- seq(tranges[i,1],tranges[i,2],length=nt)
      t2vec <- seq(tranges[j,1],tranges[j,2],length=nt)
      mtmp <- matrix(nrow=nt,ncol=nt)
      for (t1 in seq_along(t1vec)) {
        for (t2 in seq_along(t2vec)) {
          mtmp[t1,t2] <- fun(mkpar(params,c(t1vec[t1],t2vec[t2]),
                                   c(i,j)))
        }
      }
      slices[[i]][[j]] <- data.frame(var1=pn[i],var2=pn[j],
                                     expand.grid(x=t1vec,y=t2vec),
                                     z=c(mtmp))
    }
  }
  r <- list(slices=slices,ranges=tranges,params=params,dim=2)
  class(r) <- "slice"
  r
}

## flatten slice:
##  do.call(rbind,lapply(slices,do.call,what=rbind))

slices_apply <- function(s,FUN,...) {
  for (i in seq_along(s)) {
    for (j in seq_along(s[[i]])) {
      if (!is.null(s[[i]][[j]])) {
        s[[i]][[j]] <- FUN(s[[i]][[j]],...)
      }
    }
  }
  s
}

xyplot.slice <- function(x,data,type="l",scale.min=TRUE,...) {
  allslice <- do.call(rbind,x$slices)
  ## correct ordering
  allslice$var1 <- factor(allslice$var1,
                          levels=unique(as.character(allslice$var1)))
  if (scale.min) allslice$z <- allslice$z-min(allslice$z)
  pfun <- function(x1,y1,...) {
      panel.xyplot(x1,y1,...)
      if (is.null(x$params2)) {
          ## regular 1D slice
          panel.abline(v=x$params[panel.number()],col="gray")
      } else {
          ## 'transect' slice
          panel.abline(v=c(0,1),col="gray")
          panel.abline(h=y1[x1 %in% c(0,1)],col="gray")
      }
  }
  xyplot(z~x|var1,data=allslice,type=type,
         scales=list(x=list(relation="free")),
         panel=pfun,...)
}

splom.slice <- function(x,
                        data,
                        scale.min=TRUE,
                        at=NULL,
                        which.x=NULL,
                        which.y=NULL,
                        dstep=4,
                        contour=FALSE,...) {
  if (x$dim==1) stop("can't do splom on 1D slice object")
  smat <- t(x$ranges[,1:2])
  if (scale.min) {
    ## FIXME: something more elegant to flatten slice list?
    all.z <- unlist(sapply(x$slices,
                    function(x) {
                      sapply(x,
                             function(x)
                             if (is.null(x)) NULL else x[["z"]])
                    }))
    min.z <- min(all.z[is.finite(all.z)])
    ## round up to next multiple of 'dstep'
    max.z <- dstep * ((max(all.z[is.finite(all.z)])-
                       min.z) %/% dstep + 1)
    if (missing(at)) {
      at <- seq(0,max.z,by=dstep)
    }
    scale.z <- function(X) {
        X$z <- X$z-min.z
        X
    }
    x$slices <- slices_apply(x$slices,scale.z)
  }
  up0 <- function(x1, y, groups, subscripts, i, j, ...) {
      ## browser()
    sl <- x$slices[[j]][[i]]
    with(sl,panel.levelplot(x=x,y=y,z=z,contour=contour,
                            at=if (!is.null(at)) at else pretty(z),
                            subscripts=seq(nrow(sl))))
    panel.points(x$params[j],x$params[i],pch=16)
    mm <- matrix(sl$z,nrow=length(unique(sl$x)))
    ## FIXME: more robust ...
    wmin <- which(mm==min(mm),arr.ind=TRUE)
    xmin <- unique(sl$x)[wmin[1]]
    ymin <- unique(sl$y)[wmin[2]]
    panel.points(xmin,ymin,pch=1)
  }
  lp0 <- function(...) {
  }
  ## FIXME: use ?draw.colorkey to add a legend ...
  ## FIXME: make diagonal panel text smaller ???
  splom(smat,lower.panel=lp0,diag.panel=diag.panel.splom,
        upper.panel=up0,...)
}

  
## generic profiling code???
##  either need (1) optimizers with 'masks' or (2)
