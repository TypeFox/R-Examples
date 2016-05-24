# Copyright (C) 2000 Robert Gray
# distributed under the terms of the GNU public license
crr <-
# function for regression modeling of subdistribution functions
# arguments: 
#  ftime = vector of failure/censoring times
#  fstatus = vector with a unique code for each failure type and a
#     separate code for censored observations 
#  cov1 = (nobs x ncovs) matrix of fixed covariates
#  cov2 = matrix of covariates multiplied by functions of time; 
#     if used, often these covariates would also appear in cov1, 
#     to give a prop hazards effect plus a time interaction
#  tf = functions of time.  A function that takes a vector of times as
#     an argument and returns a matrix whose jth column is the value of 
#     the time function corresponding to the jth column of cov2 evaluated
#     at the input time vector.  At time tk, the
#     model includes the term cov2[,j]*tfs(tk)[,j] as a covariate.
#  cengroup = vector with different values for each group with 
#     a distinct censoring distribution (the censoring distribution
#     is estimated separately within these groups)
#  failcode = code of fstatus that denotes the failure type of interest
#  cencode = code of fstatus that denotes censored observations
#  subset = logical vector length(ftime) indicating which cases to include
#  na.action = function defining action to take for cases that have NA for
#     any of ftime, fstatus, cov1, cov2 cengroup, or subset.
#  gtol = iteration stops when a function of the gradient is < gtol.
#  maxiter = maximum # of iterations in Newton algorithm (0 computes
#     scores and var at init, but performs no iterations)
#  init = initial values of regression parameters
function(ftime,fstatus,cov1,cov2,tf,cengroup,failcode=1,cencode=0,
         subset,na.action=na.omit,gtol=1e-6,maxiter=10,init,variance=TRUE) {
  ## LS
  call <- match.call() 
  cov1.name <- deparse(substitute(cov1))
  cov1.vars <- cov2.vars <- NULL
  if(!missing(cov1)) 
    { cov1.vars <- colnames(as.matrix(cov1)) }
  cov2.name <- deparse(substitute(cov2))
  if(!missing(cov2)) 
    { cov2.vars <- colnames(as.matrix(cov2)) }
  ##
  d <- data.frame(ftime=ftime,fstatus=fstatus,
        cengroup=if (missing(cengroup)) rep(1,length(fstatus)) else cengroup)
  if (!missing(cov1)) {
    cov1 <- as.matrix(cov1)
    nc1 <- ncol(cov1)
    d <- cbind(d,cov1)
  } else {nc1 <- 0}
  if (!missing(cov2)) {
    cov2 <- as.matrix(cov2)
    nc2 <- ncol(cov2)
    d <- cbind(d,cov2)
  } else {nc2 <- 0}
  if (!missing(subset)) d <- d[subset,]
  tmp <- nrow(d)
  d <- na.action(d)
  nmis <- 0
  if (nrow(d) != tmp) {
    nmis <- tmp-nrow(d)
    cat(format(nmis),'cases omitted due to missing values\n')
  }
  d <- d[order(d$ftime),]
  ftime <- d$ftime
  cenind <- ifelse(d$fstatus==cencode,1,0)
  fstatus <- ifelse(d$fstatus==failcode,1,2*(1-cenind))
  ucg <- sort(unique.default(d$cengroup))
  cengroup <- match(d$cengroup,ucg)
  ncg <- length(ucg)
  uuu <- matrix(0,nrow=ncg,ncol=length(ftime))
  for (k in 1:ncg) {
    u <- do.call('survfit',list(formula=Surv(ftime,cenind)~1,data=
       data.frame(ftime,cenind,cengroup),subset=cengroup==k))
### note: want censring dist km at ftime-
    u <- approx(c(0,u$time,max(u$time)*(1+10*.Machine$double.eps)),c(1,u$surv,
       0),xout=ftime*(1-100*.Machine$double.eps),method='constant',f=0,rule=2)
    uuu[k,1:length(u$y)] <- u$y
#    u <- summary(u,times=sort(ftime*(1-.Machine$double.eps)))
#    uuu[k,1:length(u$surv)] <- u$surv
  }
  uft <- sort(unique(ftime[fstatus==1]))
  ndf <- length(uft)
  if (nc2 == 0) {
    cov1 <- as.matrix(d[,(1:nc1)+3])
    np <- nc1
    npt <- 0
    cov2 <- 0
    tfs <- 0
  } else if (nc1 == 0) {
    cov2 <- as.matrix(d[,(1:nc2)+3+nc1])
    npt <- np <- nc2
    cov1 <- 0
    tfs <- tf(uft)
  } else {
    cov1 <- as.matrix(d[,(1:nc1)+3])
    cov2 <- as.matrix(d[,(1:nc2)+3+nc1])
    npt <- nc2
    np <- nc1+nc2
    tfs <- tf(uft)
  }
### start of nr
  if (missing(init)) b <- rep(0,np)
  else b <- init
  stepf <- .5
  for (ll in 0:maxiter) {
    z <- .Fortran('crrfsv',as.double(ftime),as.integer(fstatus),
                  as.integer(length(ftime)),as.double(cov1),as.integer(np-npt),
                  as.integer(np),as.double(cov2),as.integer(npt),
                  as.double(tfs),as.integer(ndf),as.double(uuu),
                  as.integer(ncg),as.integer(cengroup),as.double(b),
                  double(1),double(np),double(np*np),double(np),double(np),
                  double(np*np),PACKAGE = "cmprsk")[15:17]
    if (max(abs(z[[2]])*pmax(abs(b),1)) < max(abs(z[[1]]),1)*gtol) {
      converge <- TRUE
      break
    }
    if (ll==maxiter) {
      converge <- FALSE
      break
    }
    h <- z[[3]]
    dim(h) <- c(np,np)
### better to guarantee a pd factorization, but
### matrix should be pd except in rare circumstances
    sc <- -solve(h,z[[2]])
    bn <- b+sc
    fbn <- .Fortran('crrf',as.double(ftime),as.integer(fstatus),
                  as.integer(length(ftime)),as.double(cov1),as.integer(np-npt),
                  as.integer(np),as.double(cov2),as.integer(npt),
                  as.double(tfs),as.integer(ndf),as.double(uuu),
                  as.integer(ncg),as.integer(cengroup),as.double(bn),
                  double(1),double(np),PACKAGE = "cmprsk")[[15]]
# backtracking loop
    i <- 0
    while (is.na(fbn) || fbn>z[[1]]+(1e-4)*sum(sc*z[[2]])) {
      i <- i+1
      sc <- sc*stepf
      bn <- b+sc
      fbn <- .Fortran('crrf',as.double(ftime),as.integer(fstatus),
                  as.integer(length(ftime)),as.double(cov1),as.integer(np-npt),
                  as.integer(np),as.double(cov2),as.integer(npt),
                  as.double(tfs),as.integer(ndf),as.double(uuu),
                  as.integer(ncg),as.integer(cengroup),as.double(bn),
                  double(1),double(np),PACKAGE = "cmprsk")[[15]]
      if (i>20) break
    }
    if (i>20) {
      converge <- FALSE
      break
    }
    b <- c(bn)
  }
  if (variance) {
    v <- .Fortran('crrvv',as.double(ftime),as.integer(fstatus),
                as.integer(length(ftime)),as.double(cov1),as.integer(np-npt),
                as.integer(np),as.double(cov2),as.integer(npt),
                as.double(tfs),as.integer(ndf),as.double(uuu),
                as.integer(ncg),as.integer(cengroup),as.double(b),
                double(np*np),double(np*np),double(np*np),
                double(length(ftime)*(np+1)),double(np),double(np*ncg),
                double(2*np),double(ncg*np),
                integer(ncg),double(ncg*np),double(ncg),
                PACKAGE = "cmprsk")[15:16]
    dim(v[[2]]) <- dim(v[[1]]) <- c(np,np)
    h0 <- v[[1]]
    h <- solve(v[[1]])
    v <- h %*% v[[2]] %*% t(h)
    r <- .Fortran('crrsr',as.double(ftime),as.integer(fstatus),
                as.integer(length(ftime)),as.double(cov1),as.integer(np-npt),
                as.integer(np),as.double(cov2),as.integer(npt),
                as.double(tfs),as.integer(ndf),as.double(uuu),
                as.integer(ncg),as.integer(cengroup),as.double(b),
                double(ndf*np),double(np),double(np),PACKAGE = "cmprsk")[[15]]
    r <- t(matrix(r,nrow=np))
##
  } else {
    v <- h <- h0 <- matrix(NA,np,np)
    r <- NULL
#    bj <- NULL
  }
  nobs <- length(ftime)
  b0 <- rep(0,length(b))
  fb0 <- .Fortran('crrf',as.double(ftime),as.integer(fstatus),
                  as.integer(length(ftime)),as.double(cov1),as.integer(np-npt),
                  as.integer(np),as.double(cov2),as.integer(npt),
                  as.double(tfs),as.integer(ndf),as.double(uuu),
                  as.integer(ncg),as.integer(cengroup),as.double(b0),
                  double(1),double(np),PACKAGE = "cmprsk")[[15]]
  bj <- .Fortran('crrfit',as.double(ftime),as.integer(fstatus),
                  as.integer(length(ftime)),as.double(cov1),as.integer(np-npt),
                  as.integer(np),as.double(cov2),as.integer(npt),
                  as.double(tfs),as.integer(ndf),as.double(uuu),
                  as.integer(ncg),as.integer(cengroup),as.double(b),
                  double(ndf),double(np),PACKAGE = "cmprsk")[[15]]
  if (nc1>0) {
    x1 <- paste(cov1.name, 1:nc1, sep="")
    if(is.null(cov1.vars)) cov1.vars <- x1
    else cov1.vars <- ifelse(cov1.vars=="",x1,cov1.vars)
  }
  if(nc2 > 0) { 
    x1 <- paste(cov2.name, 1:nc2, sep="")
    if (is.null(cov2.vars)) cov2.vars <- x1
    else cov2.vars <- ifelse(cov2.vars=="",x1,cov2.vars)
    x1 <- paste('tf',1:nc2,sep='')
    x2 <- colnames(tfs)
    if (!is.null(x2)) x1 <- ifelse(x2=="",x1,x2)
    cov2.vars <- paste(cov2.vars, x1, sep="*")
  }
  names(b) <- c(cov1.vars, cov2.vars)
    ##
  z <- list(coef=b,loglik=-z[[1]],score=-z[[2]],inf=h0,
            var=v,res=r,uftime=uft,bfitj=bj,
            tfs=as.matrix(tfs),converged=converge,call = call, n = nobs, 
            n.missing = nmis, loglik.null = -fb0,invinf=h)
  class(z) <- 'crr'
  z
}

"summary.crr" <- function(object, conf.int = 0.95, digits = max(options()$digits - 5, 2), ...) 
{
  beta <- object$coef
  se <- sqrt(diag(object$var))
  out <- list(call = object$call, converged = object$converged, 
              n = object$n, n.missing = object$n.missing, 
              loglik = object$loglik)
  tmp <- cbind(beta, exp(beta), se, beta/se, 
               signif(2 * (1 - pnorm(abs(beta)/se)), digits))
  dimnames(tmp) <- list(names(beta), c("coef", "exp(coef)", 
                        "se(coef)", "z", "p-value"))
  out$coef <- tmp
  if(conf.int) 
    { a <- (1 - conf.int)/2
      a <- c(a, 1 - a)
      z <- qnorm(a)
      tmp <- cbind(exp(beta), exp(-beta), 
                   exp(beta + z[1] * se), exp(beta + z[2] * se))
      dimnames(tmp) <- list(names(beta), c("exp(coef)", "exp(-coef)", 
                            paste(format(100*a, trim = TRUE, 
                                         scientific = FALSE, 
                                         digits = 3), "%", sep="")))
      out$conf.int <- tmp
  }
  df <- length(beta)
  logtest <- -2 * (object$loglik.null - object$loglik)
  out$logtest <- c(test = logtest, df = df)
  # out$rsq <- c(rsq = 1 - exp(-logtest/object$n), 
  #              maxrsq = 1 - exp(2 * object$loglik.null/object$n))
  class(out) <- "summary.crr"
  out
}

"print.summary.crr" <- function (x, digits = max(options()$digits - 4, 3), ...) 
{
    cat("Competing Risks Regression\n\n")
    if(!is.null(x$call)) 
      { cat("Call:\n")
        dput(x$call)
        cat("\n") 
      }
    if(!x$converged) 
      { cat("crr converged:", x$converged, "\n")
        return()
      }
    savedig <- options(digits = digits)
    on.exit(options(savedig))
    print(x$coef)
    cat("\n")
    print(x$conf.int)
    cat("\n")
    cat("Num. cases =", x$n)
    if(x$n.missing > 0) 
      cat(" (", x$n.missing, " cases omitted due to missing values)", sep="")
    cat("\n")
    # cat("Rsquare =", format(round(x$rsq["rsq"], 3)), "  (max possible =", 
    #    format(round(x$rsq["maxrsq"], 3)), ")\n")
    cat("Pseudo Log-likelihood =", x$loglik, "\n")
    cat("Pseudo likelihood ratio test = ", format(round(x$logtest["test"], 2)), 
        "  on ", x$logtest["df"], " df,",  "\n", sep = "")
#        "  p-value = ", format(x$logtest["pvalue"]), "\n", sep = "")
    invisible()
}

predict.crr <-
# for a crr object x, estimates subdistributions at covariate
# combinations given by rows of cov1 and cov2.  The terms in cov1
# cov2 must correspond exactly to the corresponding call to crr.
  function(object,cov1,cov2,...) {
    if (is.null(object$bfitj)) stop('predict requires variance=TRUE in crr')
    np <- length(object$coef)
    if (length(object$tfs)<=1) {
      if (length(object$coef)==length(cov1)) lhat <- cumsum(exp(sum(cov1*object$coef))*object$bfitj)
      else {
        cov1 <- as.matrix(cov1)
        lhat <- matrix(0,nrow=length(object$uftime),ncol=nrow(cov1))
        for (j in 1:nrow(cov1)) lhat[,j] <- cumsum(exp(sum(cov1[j,]*object$coef))*object$bfitj)
      }
    } else {
      if (length(object$coef)==ncol(as.matrix(object$tfs))) {
        if (length(object$coef)==length(cov2))
          lhat <- cumsum(exp(object$tfs %*% c(cov2*object$coef))*object$bfitj)
        else {
          cov2 <- as.matrix(cov2)
          lhat <- matrix(0,nrow=length(object$uftime),ncol=nrow(cov1))
          for (j in 1:nrow(cov2)) lhat[,j] <-
            cumsum(exp(object$tfs %*% c(cov2[j,]*object$coef))*object$bfitj)
        }
      } else {
        if (length(object$coef)==length(cov1)+length(cov2))
          lhat <- cumsum(exp(sum(cov1*object$coef[1:length(cov1)])+object$tfs %*%
                             c(cov2*object$coef[(np-length(cov2)+1):np]))*object$bfitj)
        else {
          cov1 <- as.matrix(cov1)
          cov2 <- as.matrix(cov2)
          lhat <- matrix(0,nrow=length(object$uftime),ncol=nrow(cov1))
          for (j in 1:nrow(cov1)) lhat[,j] <-
            cumsum(exp(sum(cov1[j,]*object$coef[1:ncol(cov1)])+object$tfs %*%
                       c(cov2[j,]*object$coef[(np-ncol(cov2)+1):np]))*object$bfitj)
        }
      }
    }
    lhat <- cbind(object$uftime,1-exp(-lhat))
    class(lhat) <- 'predict.crr'
    lhat
  }
plot.predict.crr <-
# plots estimated subdistributions from predict.crr
  function(x,lty=1:(ncol(x)-1),color=1,ylim=c(0,max(x[,-1])),xmin=0,xmax=max(x[,1]),...) {
  if (length(lty)<ncol(x)-1) lty <- rep(lty[1],ncol(x)-1)
  if (length(color)<ncol(x)-1) color <- rep(color[1],ncol(x)-1)
  if (xmax<max(x[,1])) x <- x[x[,1]<xmax,]
  times <- c(xmin,rep(x[,1],rep(2,nrow(x))),xmax)
  plot(c(xmin,xmax),ylim,type='n',...)
  for (j in 2:ncol(x)) lines(times,c(0,0,rep(x[,j],rep(2,nrow(x)))),lty=lty[j-1],col=color[j-1])
}
print.crr <-
# prints a summary of the crr fit x
  function(x,...) {
  cat('convergence: ',x$converged,'\n')
  cat('coefficients:\n')
  print(signif(x$coef,4),...)
  v <- sqrt(diag(x$var))
  cat('standard errors:\n')
  print(signif(v,4),...)
  v <- 2*(1-pnorm(abs(x$coef)/v))
  cat('two-sided p-values:\n')
  print(signif(v,2),...)
  invisible()
}

cuminc <- function(ftime,fstatus,group,strata,rho=0,cencode=0,subset,na.action=na.omit) {
# ftime=failure times, fstatus=variable which indicates the type
# of failure (and cens), group is the group variable, strata=
# strata variables for the tests (omit if none), rho is the 
# power of the weight function used in the tests, cencode is the value
# of fstatus which indicates that a time is censored (default is 0)
# subset = logical vector length(ftime) indicating which cases to include
# na.action = function defining action to take for cases that have NA for
#     any of ftime, fstatus, group, strata, or subset.
# output is a list giving the estmated cuminc
# functions (times, function values, variances) for each group, and
# a component
# values of the test statistics for comparing each cause among the
# groups.  the tests are stratified (if strata specified).  
# check lengths, and status of group and strata
  d <- data.frame(time=ftime,cause=fstatus,
    group=as.factor(if (missing(group)) rep(1,length(ftime)) else group),
    strata=as.factor(if (missing(strata)) rep(1,length(ftime)) else strata))
  if (!missing(subset)) d <- d[subset,]
  tmp <- nrow(d)
  d <- na.action(d)
  if (nrow(d) != tmp) cat(format(tmp-nrow(d)),'cases omitted due to missing values\n')
  no <- nrow(d)
  cg <- "  "
  nst <- length(levels(d$strata))
  d <- d[order(d$time),]
  ugg <- table(d$group)
  d$group <- factor(d$group,names(ugg)[ugg>0])
  ugg <- levels(d$group)
  censind <- ifelse(d$cause==cencode,0,1)
  uc <- table(d$cause[censind==1])
  if (is.factor(d$cause)) uclab <- names(uc)[uc>0]
  else uclab <- as.numeric(names(uc)[uc>0])
  nc <- length(uclab)
  ng <- length(ugg)
  if (ng>1) {
    ng1 <- ng-1
    ng2 <- ng*ng1/2
    v <- matrix(0,nrow=ng1,ncol=ng1)
    storage.mode(v) <- "double"
    vt <- double(ng2)
    s <- double(ng1)
  }
  pf <- vector("list",ng*nc)
  stat <- double(nc)
  l <- 0
  for (ii in 1:nc) {
    causeind <- ifelse(d$cause==uclab[ii],1,0)
    for (jj in 1:length(ugg)) {
      cg <- c(cg,paste(ugg[jj],uclab[ii]))
      l <- l+1
      cgind <- d$group==ugg[jj]
      ncg <- length(cgind[cgind])
      n2 <- length(unique(d$time[cgind & causeind==1]))
      n2 <- 2*n2+2
      tmp <- double(n2)
      z <- .Fortran("cinc",as.double(d$time[cgind]),as.integer(censind[cgind]),
              as.integer(causeind[cgind]),as.integer(ncg),
              x=tmp,f=tmp,v=tmp,PACKAGE = "cmprsk")
      pf[[l]] <- list(time=z$x,est=z$f,var=z$v)
    }
    if (ng>1) {
      causeind <- 2*censind-causeind
      z2 <- .Fortran("crstm",as.double(d$time),as.integer(causeind),
               as.integer(d$group),as.integer(d$strata),as.integer(no),
               as.double(rho),as.integer(nst),as.integer(ng),s,v,
               as.double(d$time),as.integer(causeind),as.integer(d$group),
               vt,s,vt,double((4+3*ng)*ng),integer(4*ng),PACKAGE = "cmprsk")
      stat[ii] <- -1
      a <- qr(z2[[10]])
      if (a$rank==ncol(a$qr)) {
        b <- diag(dim(a$qr)[1])
        stat[ii] <- z2[[9]]%*%qr.coef(a,b)%*%z2[[9]]
      }
    }
  }
  names(pf) <- cg[2:length(cg)]
  if (ng>1) {
    names(stat) <- uclab
    stat <- list(Tests=cbind(stat=stat,pv=1-pchisq(stat,ng-1),df=rep(ng-1,length(stat))))
    pf <- c(pf,stat)
  }
  attr(pf, "class") <- "cuminc"
  pf
}

print.cuminc <- function(x,ntp=4,maxtime,...) {
  if (!is.null(x$Tests)) {
    cat('Tests:\n')
    print(x$Tests)
    nc <- length(x)-1
  } else {
    nc <- length(x)
  }
  if (missing(maxtime)) {
    maxtime <- 0
    for (i in 1:nc) maxtime <- max(maxtime,x[[i]]$time)
  }
  tp <- pretty(c(0,maxtime),ntp+1)
  cat('Estimates and Variances:\n')
  print(timepoints(x,tp[-c(1,length(tp))]),...)
  invisible()
}

timepoints <- function(w,times) {
# w=list (see cuminc or km), times= times you want estimates for.
# output is a list with components est giving the estimates, var giving
# the variances,
  if (!is.null(w$Tests)) w <- w[names(w) != 'Tests']
  ng <- length(w)
  times <- sort(unique(times))
  nt <- length(times)
  storage.mode(times) <- "double"
  storage.mode(nt) <- "integer"
  ind <- matrix(0,ncol=nt,nrow=ng)
  oute <- matrix(NA,ncol=nt,nrow=ng)
  outv <- oute
  storage.mode(ind) <- "integer"
  slct <- rep(TRUE,ng)
  for (i in 1:ng) {
    if (is.null((w[[i]])$est)) { slct[i] <- FALSE} else {
      z <- .Fortran("tpoi",as.double(w[[i]][[1]]),
         as.integer(length(w[[i]][[1]])),ind[i,],times,nt,PACKAGE = "cmprsk")
      ind[i,] <- z[[3]]
      oute[i,ind[i,]>0] <- w[[i]][[2]][z[[3]]]
      if (length(w[[i]])>2) outv[i,ind[i,]>0] <- w[[i]][[3]][z[[3]]]
    }
  }
  dimnames(oute) <- list(names(w)[1:ng],as.character(times))
  dimnames(outv) <- dimnames(oute)
  list(est=oute[slct,,drop=FALSE],var=outv[slct,,drop=FALSE])
}

plot.cuminc <-  function(x,main=" ",curvlab,ylim=c(0,1),xlim,wh=2,xlab="Years",
ylab="Probability",lty=1:length(x),color=1,lwd = par('lwd'),...) {
# x is a list containing curves to be plotted. Each component of
# x is a list with the first component containing the x values
# and the second component the y values.  main = main title in the plot
# curvlab=curve labels (vector), wh=where curve labels are plotted
# 1=lower left 2=upper left 3=upper right 4=lower right
  if (!is.null(x$Tests)) x <- x[names(x) != 'Tests']
  nc <- length(x)
  if (length(lty) < nc) lty <- rep(lty[1],nc) else lty <- lty[1:nc]
  if (length(lwd) < nc) lwd <- rep(lwd[1],nc) else lwd <- lwd[1:nc]
  if (length(color) < nc) color <- rep(color[1],nc) else color <- color[1:nc]
  if (missing(curvlab)) {
    if (mode(names(x))=="NULL") {
      curvlab <- as.character(1:nc) }
    else curvlab <- names(x)[1:nc]
  }
  if (missing(xlim)) {
    xmax <- 0
    for (i in 1:nc) {
      xmax <- max(c(xmax,x[[i]][[1]]))
    }
    xlim <- c(0,xmax)
  }
  plot(x[[1]][[1]],x[[1]][[2]],type="n",ylim=ylim,xlim=xlim,
       main=main,xlab=xlab,ylab=ylab,bty="l",...)
  if (length(wh) != 2) {
      wh <- c(xlim[1],ylim[2])
  }
  u <- list(...)
  if (length(u)>0) {
    i <- pmatch(names(u),names(formals(legend)),0)
    do.call('legend',c(list(x=wh[1],y=wh[2],legend=curvlab,col=color,lty=lty,lwd=lwd,bty="n",bg=-999999),u[i>0]))
  } else {
    do.call('legend',list(x=wh[1],y=wh[2],legend=curvlab,col=color,lty=lty,lwd=lwd,bty="n",bg=-999999))
  }
#  legend(wh[1],wh[2],legend=curvlab,col=color,lty=lty,bty="n",bg=-999999,...)
  for (i in 1:nc) {
    lines(x[[i]][[1]],x[[i]][[2]],lty=lty[i],col=color[i],lwd=lwd[i],...)
  }
}

"[.cuminc" <- function(x,i,...) {
  x <- NextMethod("[")
  class(x) <- 'cuminc'
  x
}
