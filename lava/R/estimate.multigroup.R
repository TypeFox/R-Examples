###{{{ estimate.multigroup

##' @export
`estimate.multigroup` <- function(x, control=list(),
                                  estimator="gaussian",
                                  weight, weightname,
                                  weight2,
                                  id=NULL,
                                  silent=lava.options()$silent,
                                  quick=FALSE,
                                  param,
                                  cluster,
                                  ...) {
  cl <- match.call()
  optim <- list(
             iter.max=lava.options()$iter.max,
             trace=ifelse(lava.options()$debug,3,0),
             gamma=lava.options()$gamma,
             ngamma=lava.options()$ngamma,
	     backtrace=TRUE,
             gamma2=1,
             lambda=0.05,
             abs.tol=1e-9,
             epsilon=1e-10,
             delta=1e-10,
             S.tol=1e-6,
             stabil=FALSE,
             start=NULL,
             constrain=lava.options()$constrain,
             method=NULL,
             starterfun=startvalues0,
             information="E",
             meanstructure=TRUE,
             sparse=FALSE,
             lbound=1e-9,
             reindex=FALSE,
             tol=lava.options()$tol)


  if (!missing(param)) {
    oldparam <- lava.options()$param
    lava.options(param=param)
    on.exit(lava.options(param=oldparam))
  }
  if (!missing(cluster)) id <- cluster

  defopt <- lava.options()[]
  defopt <- defopt[intersect(names(defopt),names(optim))]
  optim[names(defopt)] <- defopt

  if (length(control)>0) {
    optim[names(control)] <- control
  }


  Debug("Start values...")
  if (!is.null(optim$start) & length(optim$start)==(x$npar+x$npar.mean)) {
    mystart <- optim$start
  } else {
    if (!silent) cat("Obtaining starting value...")
    if (is.null(control$starterfun) && lava.options()$param!="relative")
        optim$starterfun <- startvalues0
    mystart <- with(optim, starter.multigroup(x,meanstructure=meanstructure,starterfun=starterfun,silent=FALSE,fix=FALSE))
    if (!is.null(optim$start)) {
      pname <- names(optim$start)
      ppos <- parpos.multigroup(x,p=pname,mean=TRUE)
      if (any(!is.na(ppos)))
        mystart[ppos] <- optim$start[na.omit(match(attributes(ppos)$name,pname))]
    }
    if (!silent) cat("\n")
  }
  Debug(mystart)
  Debug("Constraints...")
  ## Setup optimization constraints
  lower <- rep(-Inf, x$npar);
  for (i in seq_len(x$ngroup)) {
    vpos <- sapply(x$parlist[[i]][variances(x$lvm[[i]],mean=FALSE)], function(y) as.numeric(substr(y,2,nchar(y))))
    if (length(vpos)>0)
    lower[vpos] <- optim$lbound
  }
  if (optim$meanstructure)
    lower <- c(rep(-Inf,x$npar.mean), lower)
  if (any(optim$constrain)) {
    if (length(optim$constrain)!=length(lower))
      constrained <- is.finite(lower)
    else
      constrained <- optim$constrain
    constrained <- which(constrained)
    lower[] <- -Inf
    optim$constrain <- TRUE
    mystart[constrained] <- log(mystart[constrained])
  }

  if (!missing(weight)) {
    if (is.character(weight)) {
      stweight <- weight
      weight <- list()
      for (i in seq_along(x$data)) {
        newweight <- as.matrix(x$data[[i]][,stweight])
        colnames(newweight) <- index(x$lvm[[i]])$endogenous[seq_len(ncol(newweight))]
        weight <- c(weight, list(newweight))
      }
    }
  } else {
    weight <- NULL
  }
  if (!missing(weight2)) {
    if (is.character(weight2)) {
      stweight2 <- weight2
      weight2 <- list()
      for (i in seq_along(x$data)) {
        newweight <- as.matrix(x$data[[i]][,stweight2,drop=FALSE])
        dropcol <- apply(newweight,2,function(x) any(is.na(x)))
        newweight <- newweight[,!dropcol,drop=FALSE]
        colnames(newweight) <- index(x$lvm[[i]])$endogenous[seq_len(ncol(newweight))]
        weight2 <- c(weight2, list(newweight))
      }
    }
  } else {
    weight2 <- NULL
  }

### Run hooks (additional lava plugins)
  myhooks <- gethook()
  newweight <- list()
  newweight2 <- list()
  newoptim <- newestimator <- NULL
  for (f in myhooks) {
    for ( i in seq_len(x$ngroup)) {
      res <- do.call(f, list(x=x$lvm[[i]],data=x$data[[i]],weight=weight[[i]],weight2=weight2[[i]],estimator=estimator,optim=optim))
      if (!is.null(res$x)) x$lvm[[i]] <- res$x
      if (!is.null(res$data)) x$data[[i]] <- res$data
      if (!is.null(res$weight)) newweight <- c(newweight,list(res$weight))
      if (!is.null(res$weight2)) newweight2 <- c(newweight2,list(res$weight2))
      if (!is.null(res$optim)) newoptim <- res$optim
      if (!is.null(res$estimator)) newestimator <- res$estimator
    }
    if (!is.null(newestimator)) estimator <- newestimator
    if (!is.null(newoptim)) optim <- newoptim
    if (!is.null(res$weight))
      if (!any(unlist(lapply(newweight,is.null)))) {
        weight <- newweight
      }
    if (!is.null(res$weight2))
      if (!any(unlist(lapply(newweight2,is.null)))) {
        weight2 <- newweight2
      }
  }


  checkestimator <- function(x,...) {
    ffname <- paste0(x,c("_objective","_gradient"),".lvm")
    exists(ffname[1])||exists(ffname[2])
  }
  if (!checkestimator(estimator)) { ## Try down/up-case version
    estimator <- tolower(estimator)
    if (!checkestimator(estimator)) {
      estimator <- toupper(estimator)
    }
  }

  Method <-  paste0(estimator, "_method", ".lvm")
  if (!exists(Method))
    Method <- "nlminb1"
  else
    Method <- get(Method)
  if (is.null(optim$method)) {
      optim$method <- Method
  }

  ## Check for random slopes
  xXx <- exogenous(x)
  Xfix <- FALSE
  Xconstrain <- FALSE
  xfix <- list()
  for (i in seq_len(x$ngroup)) {
    x0 <- x$lvm[[i]]
    data0 <- x$data[[i]]
    xfix0 <- colnames(data0)[(colnames(data0)%in%parlabels(x0,exo=TRUE))]
    xconstrain0 <- intersect(unlist(lapply(constrain(x0),function(z) attributes(z)$args)),manifest(x0))
    xfix <- c(xfix, list(xfix0))
    if (length(xfix0)>0) Xfix<-TRUE ## Yes, random slopes
    if (length(xconstrain0)>0) Xconstrain <- TRUE ## Yes, nonlinear regression
  }

  ## Non-linear parameter constraints involving observed variables? (e.g. nonlinear regression)
  constr <- c()
  XconstrStdOpt <- TRUE
  xconstrainM <- TRUE
  xconstrain <- c()
  if (Xconstrain)
  for (i in seq_len(x$ngroup)) {
    x0 <- x$lvm[[i]]
    data0 <- x$data[[i]]
    constr0 <- lapply(constrain(x0), function(z)(attributes(z)$args))
    xconstrain0 <- intersect(unlist(constr0), manifest(x0))
    xconstrain <- c(xconstrain, list(xconstrain0))
    if (length(xconstrain0)>0) {
      constrainM0 <- names(constr0)%in%unlist(x0$mean)
      for (i in seq_len(length(constr0))) {
        if (!constrainM0[i]) {
          if (xconstrain0%in%constr0[[i]]) {
            xconstrainM <- FALSE
          }
        }
      }
      if (xconstrainM & ((is.null(control$method) || optim$method=="nlminb0") & (lava.options()$test & estimator=="gaussian")) ) {
        XconstrStdOpt <- FALSE
        optim$method <- "nlminb0"
        if (is.null(control$constrain)) control$constrain <- TRUE
      }
    }
  }

  ## Define objective function and first and second derivatives
  ObjectiveFun  <- paste0(estimator, "_objective", ".lvm")
  GradFun  <- paste0(estimator, "_gradient", ".lvm")
  if (!exists(ObjectiveFun) & !exists(GradFun)) stop("Unknown estimator.")

  InformationFun <- paste0(estimator, "_hessian", ".lvm")

  parord <- modelPar(x,seq_len(with(x,npar+npar.mean)))$p
  mymodel <- x

  parkeep <- c()
  myclass <- c("multigroupfit","lvmfit")
  myfix <- list()

  if (Xfix |  (Xconstrain & XconstrStdOpt | !lava.options()$test)) { ## Model with random slopes:
#############################################################

    if (Xfix) {
      myclass <- c(myclass,"lvmfit.randomslope")
      for (k in seq_len(x$ngroup)) {
        x1 <- x0 <- x$lvm[[k]]
        data0 <- x$data[[k]]

        nrow <- length(vars(x0))
        xpos <- lapply(xfix[[k]],function(y) which(regfix(x0)$labels==y))
        colpos <- lapply(xpos, function(y) ceiling(y/nrow))
        rowpos <- lapply(xpos, function(y) (y-1)%%nrow+1)
        myfix0 <- list(var=xfix[[k]], col=colpos, row=rowpos)
        myfix <- c(myfix, list(myfix0))

        for (i in seq_along(myfix0$var))
          for (j in seq_along(myfix0$col[[i]]))
            regfix(x0,
                   from=vars(x0)[myfix0$row[[i]][j]],to=vars(x0)[myfix0$col[[i]][j]]) <-
                     colMeans(data0[,myfix0$var[[i]],drop=FALSE],na.rm=TRUE)
        index(x0) <- reindex(x0,zeroones=TRUE,deriv=TRUE)
        x$lvm[[k]] <- x0
        yvars <- endogenous(x0)
        parkeep <- c(parkeep, parord[[k]][coef(x1,mean=TRUE,fix=FALSE)%in%coef(x0,mean=TRUE,fix=FALSE)])
      }
      parkeep <- sort(unique(parkeep))
      ## Alter start-values:

      if (length(mystart)!=length(parkeep))
        mystart <- mystart[parkeep]
      lower <- lower[parkeep]
      x <- multigroup(x$lvm,x$data,fix=FALSE,exo.fix=FALSE)
    }

    parord <- modelPar(x,seq_along(mystart))$p
    mydata <- list()
    for (i in seq_len(x$ngroup)) {
      mydata <- c(mydata, list(as.matrix(x$data[[i]][,manifest(x$lvm[[i]])])))
    }

    myObj <- function(theta) {
      if (optim$constrain)
        theta[constrained] <- exp(theta[constrained])
      pp <- modelPar(x,theta)$p
      res <- 0
    for (k in seq_len(x$ngroup)) {
        x0 <- x$lvm[[k]]
        data0 <- x$data[[k]]
        if (Xfix) {
          xfix0 <- xfix[[k]]
          myfix0 <- myfix[[k]]
        }
        p0 <- pp[[k]]
        myfun <- function(ii) {
          if (Xfix)
          for (i in seq_along(myfix0$var)) {
            x0$fix[cbind(myfix0$row[[i]],myfix0$col[[i]])] <-
              index(x0)$A[cbind(myfix0$row[[i]],myfix0$col[[i]])] <-
                data0[ii,xfix0[i]]
          }
          if (is.list(weight2[[k]][ii,])) {
            res <- do.call(ObjectiveFun, list(x=x0, p=p0,
                                              data=data0[ii,manifest(x0),drop=FALSE],
                                              n=1, S=NULL, weight=weight[[k]][ii,],
                                              weight2=weight2[[k]]))

          } else {
            res <- do.call(ObjectiveFun, list(x=x0, p=p0,
                                              data=data0[ii,manifest(x0),drop=FALSE],
                                              n=1, S=NULL, weight=weight[[k]][ii,],
                                              weight2=weight2[[k]][ii,]))
          }
          return(res)
        }
        res <- res + sum(sapply(seq_len(nrow(mydata[[k]])),myfun))
      }
      res
    }

    myGrad <- function(theta) {
      if (optim$constrain) {
        theta[constrained] <- exp(theta[constrained])
      }
      pp <- modelPar(x,theta)$p
      D0 <- res <- rbind(numeric(length(mystart)))
      for (k in seq_len(x$ngroup)) {
        if (Xfix) {
          myfix0 <- myfix[[k]]
        }
        x0 <- x$lvm[[k]]
        myfun <- function(ii) {
          if (Xfix)
          for (i in seq_along(myfix0$var)) {
            x0$fix[cbind(myfix0$row[[i]],myfix0$col[[i]])] <-
              index(x0)$A[cbind(myfix0$row[[i]],myfix0$col[[i]])] <-
                x$data[[k]][ii,xfix[[k]][i]]
          }
          if (is.list(weight2[[k]][ii,])) {

          } else {
            val <- do.call(GradFun, list(x=x0, p=pp[[k]],
                                         data=mydata[[k]][ii,,drop=FALSE], n=1,
                                         S=NULL,
                                         weight=weight[[k]][ii,],
                                         weight2=weight2[[k]][ii,]))
          }
          return(val)
        }
        D <- D0; D[parord[[k]]] <- rowSums(sapply(seq_len(nrow(mydata[[k]])),myfun))
        res <- res+D
      }
      if (optim$constrain) {
        res[constrained] <- res[constrained]*theta[constrained]
      }
      return(as.vector(res))
    }

    myInformation <- function(theta) {
      theta0 <- theta
      if (optim$constrain) {
        theta[constrained] <- exp(theta[constrained])
      }
      pp <- modelPar(x,theta)$p
      I0 <- res <- matrix(0,length(theta),length(theta))
      grad <- grad0 <- numeric(length(theta))
      for (k in seq_len(x$ngroup)) {
        x0 <- x$lvm[[k]]
        if (Xfix) {
          myfix0 <- myfix[[k]]
        }
        myfun <- function(ii) {
          if (Xfix)
          for (i in seq_along(myfix0$var)) {
            x0$fix[cbind(myfix0$row[[i]],myfix0$col[[i]])] <- index(x0)$A[cbind(myfix0$row[[i]],myfix0$col[[i]])] <-
              x$data[[k]][ii,xfix[[k]][i]]
          }
          I <- I0
          J <- do.call(InformationFun,
                       list(x=x0, p=pp[[k]],
                            data=mydata[[k]][ii,], n=1,
                            S=NULL,
                            weight=weight[[k]][ii,],
                            weight2=weight2[[k]][ii,],
                            type=optim$information
                            )
                       )
          D <- grad0
          if (!is.null(attributes(J)$grad)) {
            D[ parord[[k]] ] <- attributes(J)$grad
            attributes(I)$grad <- D
          }
          I[ parord[[k]], parord[[k]] ] <- J
          return(I)
        }
        L <- lapply(seq_len(nrow(x$data[[k]])),function(x) myfun(x))
        if (!is.null(attributes(L[[1]])$grad))
          grad <- grad + rowSums(matrix((unlist(lapply(L,function(x) attributes(x)$grad))),ncol=length(L)))
        res <- res + apply(array(unlist(L),dim=c(length(theta),length(theta),nrow(x$data[[k]]))),c(1,2),sum)
      }
      if (!is.null(attributes(L[[1]])$grad))
        attributes(res)$grad <- grad
      return(res)
    }
  } else { ## Model without random slopes:
###########################################################


    ## Non-linear parameter constraints involving observed variables? (e.g. nonlinear regression)
    yconstrain <- c()
    iconstrain <- c()
    xconstrain <- c()
    for (j in seq_len(x$ngroup)) {
      x0 <- x$lvm[[j]]
      data0 <- x$data[[j]]
      xconstrain0 <- c()
      for (i in seq_len(length(constrain(x0)))) {
        z <- constrain(x0)[[i]]
        xx <- intersect(attributes(z)$args,manifest(x0))
        if (length(xx)>0) {
          warg <- setdiff(attributes(z)$args,xx)
          wargidx <- which(attributes(z)$args%in%warg)
          exoidx <- which(attributes(z)$args%in%xx)
          parname <- names(constrain(x0))[i]
          y <- names(which(unlist(lapply(intercept(x0),function(x) x==parname))))
          el <- list(i,y,parname,xx,exoidx,warg,wargidx,z)
          names(el) <- c("idx","endo","parname","exo","exoidx","warg","wargidx","func")
          xconstrain0 <- c(xconstrain0,list(el))
        }
      }
      yconstrain0 <- unlist(lapply(xconstrain0,function(x) x$endo))
      iconstrain0 <- unlist(lapply(xconstrain0,function(x) x$idx))
      xconstrain <- c(xconstrain, list(xconstrain0))
      yconstrain <- c(yconstrain, list(yconstrain0))
      iconstrain <- c(iconstrain, list(iconstrain0))
    }

    MkOffset <- function(pp,x,data,xconstrain,grad=FALSE) {
      if (length(xconstrain)>0) {
        Mu <- matrix(0,nrow(data),length(vars(x))); colnames(Mu) <- vars(x)
        M <- modelVar(x,p=pp,data=data)
        M$parval <- c(M$parval,  x$mean[unlist(lapply(x$mean,is.numeric))])
        for (i in seq_len(length(xconstrain))) {
          pp <- unlist(M$parval[xconstrain[[i]]$warg]);
          myidx <- with(xconstrain[[i]],order(c(wargidx,exoidx)))
          mu <- with(xconstrain[[i]],
                     apply(data[,exo,drop=FALSE],1,
                           function(x) func(
                                         unlist(c(pp,x))[myidx])))
          Mu[,xconstrain[[i]]$endo] <- mu
        }
        offsets <- Mu%*%t(M$IAi)[,endogenous(x)]
        return(offsets)
      }
      return(NULL)
    }


    myObj <- function(theta) {
      theta0 <- theta
      if (optim$constrain) {
        theta[constrained] <- exp(theta[constrained])
      }
      pp <- modelPar(x,theta)$p
      res <- c()
      for (i in seq_len(x$ngroup)) {
        offset <- MkOffset(pp[[i]],x$lvm[[i]],x$data[[i]],xconstrain[[i]])
        x0 <- x$lvm[[i]]
        data0 <- x$data[[i]][,index(x$lvm[[i]])$manifest,drop=FALSE]
        S <- x$samplestat[[i]]$S
        mu <- x$samplestat[[i]]$mu
        n <- x$samplestat[[i]]$n
        if (!is.null(offset)) {
          x0$constrain[iconstrain[[i]]] <- NULL
          pd <- procdata.lvm(x0,data0[,endogenous(x0),drop=FALSE]-offset)
          S[endogenous(x0),endogenous(x0)] <- pd$S
          mu[endogenous(x0)] <- pd$mu
          n <- pd$n
          x0$mean[yconstrain[[i]]] <- 0
        }
        res <- c(res,
                 do.call(ObjectiveFun, list(x=x0, p=pp[[i]], data=data0, S=S, mu=mu, n=n, weight=weight[[i]], weight2=weight2[[i]], offset=offset)))

      }
        sum(res)
    }

    if (!exists(GradFun)) {
      myGrad <- NULL
    } else  {
      myGrad <- function(theta) {
        theta0 <- theta
        if (optim$constrain) {
          theta[constrained] <- exp(theta[constrained])
        }
        pp <- modelPar(x,theta)$p
        D0 <- res <- rbind(numeric(length(theta)))
        for (i in seq_len(x$ngroup)) {
          repval <- with(x$samplestat[[i]],
                         do.call(GradFun, list(x=x$lvm[[i]],p=pp[[i]],
                                               data=x$data[[i]][,index(x$lvm[[i]])$manifest,drop=FALSE],
                                               S=S,mu=mu,n=n,
                                               weight=weight[[i]], weight2=weight2[[i]])))
          D <- D0; D[ parord[[i]] ] <- repval
        res <- res + D
        }
        if (optim$constrain) {
          res[constrained] <- res[constrained]*theta[constrained]
        }
        return(as.vector(res))
      }
    }

    myInformation <- function(theta) {
      theta0 <- theta
      if (optim$constrain) {
        theta[constrained] <- exp(theta[constrained])
      }
      pp <- modelPar(x,theta)$p
      I0 <- res <- matrix(0,length(theta),length(theta))
      for (i in seq_len(x$ngroup)) {
        I <- I0;
        I[ parord[[i]], parord[[i]] ] <- with(x$samplestat[[i]], do.call(InformationFun, list(p=pp[[i]], x=x$lvm[[i]], data=x$data[[i]],
                                                                                              S=S, mu=mu, n=n, weight=weight[[i]],
                                                                                              weight2=weight2[[i]],
                                                                                              type=optim$information)))
        res <- res + I
      }
      D <- myGrad(theta0)
      if (optim$constrain) {
        res[constrained,-constrained] <- apply(res[constrained,-constrained,drop=FALSE],2,function(x) x*theta[constrained]);
        res[-constrained,constrained] <- t(res[constrained,-constrained])
        if (sum(constrained)==1) {
          res[constrained,constrained] <- res[constrained,constrained]*outer(theta[constrained],theta[constrained]) - (D[constrained])
        } else {
          res[constrained,constrained] <- res[constrained,constrained]*outer(theta[constrained],theta[constrained]) - diag(D[constrained],nrow=length(constrained))
        }
      }
      attributes(res)$grad <- D
      return(res)
    }
  }

##############################################################


  if (!exists(InformationFun)) myInformation <- NULL
  else if (is.null(get(InformationFun))) myInformation <- NULL
  if (is.null(get(GradFun))) myGrad <- NULL

  if (!silent) cat("Optimizing objective function...\n")
  if (lava.options()$debug) {
    print(lower)
    print(optim$constrain)
    print(optim$method)
  }
  opt <- do.call(optim$method,
                 list(start=mystart, objective=myObj, gradient=myGrad, hessian=myInformation, lower=lower, control=optim))
##  if (!silent) cat("\n")

  opt$estimate <- opt$par
  if (optim$constrain) {
    opt$estimate[constrained] <- exp(opt$estimate[constrained])
  }
  if (quick) return(list(opt=opt,vcov=NA))

  if (is.null(myGrad) | !XconstrStdOpt ) {
    ## if (!requireNamespace("numDeriv")) {
    ##   opt$gradient <- naiveGrad(myObj, opt$estimate)
    ## } else {
      opt$gradient <- numDeriv::grad(myObj, opt$par, method=lava.options()$Dmethod)
  } else {
      opt$gradient <- myGrad(opt$estimate)
  }

  if (!is.null(opt$convergence)) {
      if (opt$convergence!=0) warning("Lack of convergence. Increase number of iteration or change starting values.")
  } else if (!is.null(opt$gradient) && mean(opt$gradient)^2>1e-3) warning("Lack of convergence. Increase number of iteration or change starting values.")

  if (!XconstrStdOpt) {
    myInformation <- function(theta) information(x,p=theta)
  } else {
  if (is.null(myInformation)) {
##     if (!requireNamespace("numDeriv")) stop("I do not know how to calculate the asymptotic variance of this estimator.
## For numerical approximation please install the library 'numDeriv'.")
    if (!is.null(myGrad) & XconstrStdOpt)
      myInformation <- function(theta) numDeriv::jacobian(myGrad, theta, method=lava.options()$Dmethod)
    else {
      myInformation <- function(theta) numDeriv::hessian(myObj, theta)
    }
  }
}
  I <- myInformation(opt$estimate)
  asVar <- tryCatch(Inverse(I),
                    error=function(e) matrix(NA, length(mystart), length(mystart)))
    
  res <- list(model=x, model0=mymodel, call=cl, opt=opt, meanstructure=optim$meanstructure, vcov=asVar, estimator=estimator, weight=weight, weight2=weight2, cluster=id)
  class(res) <- myclass

  myhooks <- gethook("post.hooks")
  for (f in myhooks) {
    res0 <- do.call(f,list(x=res))
    if (!is.null(res0))
      res <- res0
  }

  return(res)
}

###}}}

###{{{ estimate.list

estimate.lvmlist <-
function(x, data, silent=lava.options()$silent, fix, missing=FALSE,  ...) {
  if (base::missing(data)) {
    return(estimate(x[[1]],x[[2]],missing=missing,...))
  }
  nm <- length(x)
  if (nm==1) {
    return(estimate(x[[1]],data,missing=missing,...))
  }
  if (!all(unlist(lapply(x, function(y) inherits(y,"lvm"))))) stop ("Expected a list of 'lvm' objects.")
  if (is.data.frame(data)) {
    warning("Only one dataset - going for standard analysis on each submodel.")
    res <- c()
    for (i in seq_len(nm)) {
      res <- c(res, list(estimate(x[[i]],data=data,silent=TRUE,missing=missing, ...)))
    }
    return(res)
  }

  if (nm!=length(data)) stop("Supply dataset for each model")

  Xfix <- FALSE
  xfix <- list()
  for (i in seq_along(x)) {
    data0 <- data[[i]]
    xfix0 <- colnames(data0)[(colnames(data0)%in%parlabels(x[[i]],exo=TRUE))]
    xfix <- c(xfix, list(xfix0))
    if (length(xfix0)>0) { ## Yes, random slopes
      Xfix<-TRUE
    }
  }
  if (base::missing(fix)) {
    fix <- ifelse(Xfix,FALSE,TRUE)
  }


  mg <- multigroup(x,data,fix=fix,missing=missing,...)
  res <- estimate(mg,...)

  return(res)
}

###}}}
