
get.emdbook.packages <- function() {
   pkglist = c(## "adapt",
     "bbmle","chron",
     "coda","ellipse","ggplot2",
     "gplots","gtools","gdata",
     "MCMCpack","deSolve","plotrix","R2WinBUGS","reshape","rgl",
     "scatterplot3d")
   inst.pkgs = rownames(installed.packages())
   newpkgs <- pkglist[!pkglist %in% inst.pkgs]
   if (length(newpkgs)>0) {
      install.packages(newpkgs)
   ##  do.call("install.packages",list(pkglist))
   ## sapply(pkglist,install.packages)
   }
   warning("The adapt package is no longer available. You can work through 99% of the material in _Ecological Models and Data_ without it; for more information see http://emdbolker.wikidot.com/r")
 }


## expr: expression to evaluate (raw form)
## meanval: values of the mean (possibly named)
## vars: names
## Sigma: var-cov function
deltavar <- function(fun,meanval=NULL,vars,Sigma,verbose=FALSE) {
  expr <- as.expression(substitute(fun))
  nvals <- length(eval(expr,envir=as.list(meanval)))
  vecexp <- nvals>1 ## is the result a vector?
  if (missing(vars)) {
    if (missing(meanval) || is.null(names(meanval)))
      stop("must specify either variable names or named values for means")
    vars <- names(meanval)
  }
  derivs <- try(lapply(vars,D,expr=expr),silent=TRUE)
  symbderivs <- TRUE
  if (inherits(derivs,"try-error")) {
      if (length(grep("is not in the derivatives table",derivs))) {
          ## take numeric derivative
          symbderivs <- FALSE
          warning("some symbols not in derivative table, using numeric derivatives")
          nderivs <- with(as.list(meanval),
                         numericDeriv(expr[[1]],theta=vars))
          nderivs <- attr(nderivs,"gradient")
      } else {
          stop(paste("Error within derivs:",derivs))
      }
  } else {
      nderivs <- sapply(derivs,eval,envir=as.list(meanval))
  }
  if (verbose) {
      if (symbderivs) {
          cat("symbolic derivs:\n")
          print(derivs)
      }
      cat("value of derivs:\n")
      print(nderivs)
  }
  if (!is.matrix(Sigma) && length(Sigma)>1) Sigma <- diag(Sigma)
  ## if (!is.matrix(Sigma)) sum(Sigma*nderivs^2) else
  if (vecexp && is.list(nderivs)) nderivs <- do.call("cbind",nderivs)
  if (is.matrix(nderivs)) {
    r <- apply(nderivs,1,function(z) c(z %*% Sigma %*% matrix(z)))
  } else r <- c(nderivs %*% Sigma %*% matrix(nderivs))
  ## really only want diagonal
  r
}

deltamethod <- function(fun,z,var="x",params=NULL,max.order=2) {
  d0 <- as.expression(substitute(fun))
  dvals <- list()
  dvals[[1]] <- d0
  for (i in 2:(max.order+1)) {
    dvals[[i]] <- D(dvals[[i-1]],var)
  }
  mvals <- numeric(max.order-1)
  m <- mean(z)
  for (i in 1:(max.order-1)) {
    mvals[i] <- mean((z-m)^(i+1))
  }
  mvals[1] <- var(z)  ## kluge
  ev1 <- c(as.list(params),list(m))
  ev2 <- c(as.list(params),list(z))
  names(ev1)[length(ev1)] <-   names(ev2)[length(ev2)] <- var
  r0 <- mean(eval(d0,ev2)) ## true value
  evals <- sapply(dvals,eval,ev1)  ## evaluated derivatives
  gvals <- gamma(c(1,3:(max.order+1)))
  deltavals <- cumsum(c(1,mvals)*c(evals[-(2)])/gvals)
  results <- c(r0,deltavals)
  names(results) = c("delta","E(f(x))",paste("delta",2:max.order,sep=""))
  results
}

lseq <- function(from,to,length.out) {
  exp(seq(log(from),log(to),length.out=length.out))
}

## utility function for formatting
scinot <- function(x,format=c("latex","expression"),delim="$",
                   pref="",...) {
  format <- match.arg(format)
  y <- strsplit(as.character(formatC(x,format="e",...)),"e")[[1]]
  y[1] <- gsub("^0+","",y[1])
  y[2] <- ifelse(length(grep("^\\+",y[2]))>0,
                 gsub("^\\+0+","",y[2]),
                 gsub("^-0+","-",y[2]))
  if (format=="latex") { 
    v <- paste(delim,y[1],"\\\\times 10^{",y[2],"}",delim,sep="")
  } else if (format=="expression") {
    if (as.numeric(y[1])==1) {
      v <- substitute(expression(paste(pref,10^b)),list(pref=pref,b=as.numeric(y[2])))
    } else {
      v <- substitute(expression(paste(pref,a %*% 10^b)),
                      list(pref=pref,a=as.numeric(y[1]),b=as.numeric(y[2])))
    }
  }
  v
}

## axis in scientific notation/expressions
## obsolete via sfsmisc?
axis.scinot <- function(side,at) {
  if (missing(at)) at <- axTicks(side)
  axis(side=side,labels=FALSE)
  invisible(lapply(at,
         function(a) mtext(side=side,at=a,eval(scinot(a,"expression")),
                           line=par("mgp")[2])))
}

## convert R2WinBUGS output to coda/mcmc
as.mcmc.bugs <- function(x) {
  if (x$n.chains>1) {
    z <- list()
    for (i in 1:x$n.chains) {
      z[[i]] <- mcmc(x$sims.array[,i,],start=1,thin=x$n.thin)
    }
    class(z) <- "mcmc.list"
  } else {
    z <- mcmc(x$sims.matrix,start=1,thin=x$n.thin)
  } 
  return(z)
}

## credible interval for a theoretical distribution
tcredint <- function(dist,parlist,ranges,level=0.95,eps=1e-5,verbose=FALSE) {
  qfun = function(x) do.call(paste("q",dist,sep=""),c(list(x),parlist))
  dfun = function(x) do.call(paste("d",dist,sep=""),c(list(x),parlist))
  pfun = function(x) do.call(paste("p",dist,sep=""),c(list(x),parlist))
  if (missing(ranges))  ## set upper/lower limits for search by quantiles
    ranges <- qfun(c(eps,0.5,1-eps))
  lims <- function(pdens) { ## find lower and upper values for which prob dens = target value
    lower <- uniroot(function(p) {dfun(p)-pdens},
                     interval=ranges[1:2])$root
    upper <- uniroot(function(p) {dfun(p)-pdens},
                     interval=ranges[2:3])$root
    c(lower,upper)
  }
  limarea <- function(pdens) { ## find area between target values
     intlim <- lims(pdens)
     d <- diff(pfun(intlim))
     ## cat(pdens,intlim,d,"\n")
     d
   }
  ## these limits must be within limits set above
  v1 <- qfun(c(0.6,0.9999)) # quantiles
  v2 <- dfun(v1)  ## bracketing densities
  u <- uniroot(function(x) {limarea(x)-level},
                 interval=v2)
  intlim <- lims(u$root)
  r = c(intlim,dfun(intlim[1]),limarea(u$root))
  names(r) = c("lower","upper","p","area")
  if (verbose) r else r[1:2]
}

## credible interval for (1D) posterior distribution stored in a numeric
## vector.  assumed unimodal!  pvec (vector of parameter values),
## npost (vector of posterior densities), level, tolerance
ncredint <- function(pvec,npost,level=0.95,tol=0.01,verbose=FALSE) {
  dx = diff(pvec)[1]
  cumdist <- cumsum(npost)*dx
  midpt <- which.min(abs(cumdist-0.5))
  lims <- function(pdens) { ## find lower and upper values for which
    ## prob dens is closest to target value
    lower <- which.min(abs(npost[1:midpt]-pdens))
    upper <- which.min(abs(npost[(midpt+1):length(npost)]-pdens))+midpt
    c(lower,upper)
  }
  limarea <- function(pdens) {
    intlim <- lims(pdens)
    d <- sum(npost[intlim[1]:intlim[2]])*dx
    ##    cat(pdens,intlim,d,"\n")
    d
  }
  ## find credible interval
  v2 <- seq(0,max(npost),by=tol)
  vals <- sapply(v2,limarea)
  w <- which.min(abs(vals-level))
  r = c(pvec[lims(v2[w])],v2[w],limarea(v2[w]))
  names(r) = c("lower","upper","p","area")
  ## credible intervals; posterior density, area
  if (verbose) return(r) else return(r[1:2])
}



calcslice <- function(fit1,fit2,fn=fit1@minuslogl,
                      range=c(-0.1,1.1),
                      n=400) {
  slicep = seq(range[1],range[2],length=n)
  slicepars = t(sapply(slicep,function(x) (1-x)*coef(fit1)+x*coef(fit2)))
  ## FIXME: warning about parnames from R CMD check
  if (is.null(parnames(fn))) {
    dd <- fit1@data
    ff <- names(formals(fn))
    if (!("..." %in% ff)) {
      dd <- dd[names(dd) %in% ff]
    }
    v = apply(slicepars,1,function(x) do.call(fn,c(as.list(x,dd))))
  } else { ## vector-argument function
    v = apply(slicepars,1,fn)
  }
  list(x=slicep,y=v)
}


trcoef <- function(x,inverse=FALSE) {
  ## n.b. logit should come before log
  vals=list(pattern=c("logit","log","sqrt"),
    fun=c(plogis,exp,function(x) {x^2}),
    ifun=c(qlogis,log,sqrt))
  if (is.null(attr(x,"transf"))) attr(x,"transf") = character(length(x))
  for (p in 1:length(vals[[1]])) {
    pat = vals$pattern[[p]]
    fun = vals$fun[[p]]
    ifun = vals$ifun[[p]]
    if (!inverse) {
      w = grep(paste("^",pat,sep=""),names(x))
      if (length(w)>0) {
        attr(x,"transf")[w]=pat
        names(x)[w] = gsub(pat,"",names(x)[w])      
        x[w] = sapply(x[w],fun)
      }
    } else {
      w = which(attr(x,"transf")==pat)
      if (length(w)>0) {
        attr(x,"transf")[w] = ""
        names(x)[w] = paste(pat,names(x)[w],sep="")
        x[w] = sapply(x[w],ifun)
      }
    }
  }
  if (all(attr(x,"transf")=="")) attr(x,"transf") = NULL
  x
}

