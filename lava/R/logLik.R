###{{{ logLik.lvm

##' @export
logLik.lvm <- function(object,p,data,model="gaussian",indiv=FALSE,S,mu,n,debug=FALSE,weight=NULL,weight2=NULL,...) {
  cl <- match.call()
  xfix <- colnames(data)[(colnames(data)%in%parlabels(object,exo=TRUE))]

  constr <- lapply(constrain(object), function(z)(attributes(z)$args))
  xconstrain <- intersect(unlist(constr), manifest(object))
  xconstrainM <- TRUE
  if (length(xconstrain)>0) {
    constrainM <- names(constr)%in%unlist(object$mean)
    for (i in seq_len(length(constr))) {
      if (!constrainM[i]) {
        if (constr[[i]]%in%xconstrain) xconstrainM <- FALSE
      }
    }
  }

  Debug(xfix,debug)
  if (missing(n)) {
    n <- nrow(data)
    if (is.null(n)) n <- data$n
  }
  lname <- paste0(model,"_logLik.lvm")
  logLikFun <- get(lname)


  if (length(xfix)>0 | (length(xconstrain)>0 & !xconstrainM & !lava.options()$test & model!="gaussian")) { ##### Random slopes!
    x0 <- object
    if (length(xfix)>0) {
      Debug("random slopes...",debug)
      nrow <- length(vars(object))
      xpos <- lapply(xfix,function(y) which(regfix(object)$labels==y))
      colpos <- lapply(xpos, function(y) ceiling(y/nrow))
      rowpos <- lapply(xpos, function(y) (y-1)%%nrow+1)
      myfix <- list(var=xfix, col=colpos, row=rowpos)
      for (i in seq_along(myfix$var))
        for (j in seq_along(myfix$col[[i]])) {
          regfix(x0, from=vars(x0)[myfix$row[[i]][j]],to=vars(x0)[myfix$col[[i]][j]]) <-
            data[1,myfix$var[[i]]]
        }
      index(x0) <- reindex(x0,zeroones=TRUE,deriv=TRUE)
    }
    k <- length(index(x0)$manifest)
    myfun <- function(ii) {
      if (length(xfix)>0)
        for (i in seq_along(myfix$var)) {
          index(x0)$A[cbind(myfix$row[[i]],myfix$col[[i]])] <- data[ii,myfix$var[[i]]]
        }
      return(logLikFun(x0,data=data[ii,,drop=FALSE], p=p,weight=weight[ii,,drop=FALSE],weight2=weight2[ii,,drop=FALSE],
                       model=model,debug=debug,indiv=indiv,...))
    }
    loglik <- sapply(seq_len(nrow(data)),myfun)
    if (!indiv) {
      loglik <- sum(loglik)
      n <- nrow(data)
      attr(loglik, "nall") <- n
      attr(loglik, "nobs") <- n
      attr(loglik, "df") <- length(p)
      class(loglik) <- "logLik"
    }
    return(loglik)
  }

  if (xconstrainM) {
    xconstrain <- c()
    for (i in seq_len(length(constrain(object)))) {
      z <- constrain(object)[[i]]
      xx <- intersect(attributes(z)$args,manifest(object))
      if (length(xx)>0) {
        warg <- setdiff(attributes(z)$args,xx)
        wargidx <- which(attributes(z)$args%in%warg)
        exoidx <- which(attributes(z)$args%in%xx)
        parname <- names(constrain(object))[i]
        y <- names(which(unlist(lapply(intercept(object),function(x) x==parname))))
        el <- list(i,y,parname,xx,exoidx,warg,wargidx,z)
        names(el) <- c("idx","endo","parname","exo","exoidx","warg","wargidx","func")
        xconstrain <- c(xconstrain,list(el))
      }
    }
    if (length(xconstrain)>0) {
      yconstrain <- unlist(lapply(xconstrain,function(x) x$endo))
      iconstrain <- unlist(lapply(xconstrain,function(x) x$idx))

      Mu <- matrix(0,nrow(data),length(vars(object))); colnames(Mu) <- vars(object)
      M <- modelVar(object,p=p,data=data)
      M$parval <- c(M$parval,  object$mean[unlist(lapply(object$mean,is.numeric))])
      for (i in seq_len(length(xconstrain))) {
        pp <- unlist(M$parval[xconstrain[[i]]$warg]);
        myidx <- with(xconstrain[[i]],order(c(wargidx,exoidx)))
        mu <- with(xconstrain[[i]],
                   apply(data[,exo,drop=FALSE],1,
                         function(x) {
                          func(unlist(c(pp,x))[myidx])
                        }))
        Mu[,xconstrain[[i]]$endo] <- mu
      }
      offsets <- Mu%*%t(M$IAi)[,endogenous(object),drop=FALSE]
      object$constrain[iconstrain] <- NULL
      object$mean[yconstrain] <- 0
      loglik <- do.call(lname, c(list(object=object,p=p,data=data,indiv=indiv,weight=weight,weight2=weight2,offset=offsets),list(...)))
    } else {
      cl[[1]] <- logLikFun
      loglik <- eval.parent(cl)
    }
  } else {
    cl[[1]] <- logLikFun
    loglik <- eval.parent(cl)
  }

  if (is.null(attr(loglik,"nall")))
    attr(loglik, "nall") <- n
  if (is.null(attr(loglik,"nobs")))
    attr(loglik, "nobs") <- n##-length(p)
  if (is.null(attr(loglik,"df")))
    attr(loglik, "df") <- length(p)
  class(loglik) <- "logLik"
  return(loglik)
}

###}}}

###{{{ gaussian_loglik

##' @export
gaussian_logLik.lvm <- function(object,p,data,
                          type=c("cond","sim","exo","sat","cond2"),
                          weight=NULL, indiv=FALSE, S, mu, n, offset=NULL, debug=FALSE, meanstructure=TRUE,...) {
  exo.idx <- with(index(object), exo.obsidx)
  endo.idx <- with(index(object), endo.obsidx)
  if (type[1]=="exo") {
    if (length(exo.idx)==0 || is.na(exo.idx))
      return(0)
  }

  cl <- match.call()
  if (type[1]=="cond") {
    cl$type <- "sim"
    L0 <- eval.parent(cl)
    cl$type <- "exo"
    L1 <- eval.parent(cl)
    loglik <- L0-L1
    return(loglik)
  }

  if (missing(n)) {
    if (is.vector(data)) n <- 1
    else n <- nrow(data)
  }
  k <- length(index(object)$manifest)

  if (!is.null(offset) && type[1]!="exo") {
    data[,colnames(offset)] <- data[,colnames(offset)]-offset
  }

  if (type[1]=="sat") {
    if (missing(S)) {
      d0 <- procdata.lvm(object,data=data)
      S <- d0$S; mu <- d0$mu; n <- d0$n

    }
    if (missing(p)) p <- rep(1,length(coef(object)))
    L1 <- logLik(object,p,data,type="exo",meanstructure=meanstructure)
    ##    Sigma <- (n-1)/n*S ## ML = 1/n * sum((xi-Ex)^2)
    Sigma <- S
    loglik <- -(n*k)/2*log(2*base::pi) -n/2*(log(det(Sigma)) + k) - L1
    P <- length(endo.idx)
    k <- length(exo.idx)
    npar <- P*(1+(P-1)/2)
    if (meanstructure) npar <- npar+ (P*k + P)
    attr(loglik, "nall") <- n
    attr(loglik, "nobs") <- n
    attr(loglik, "df") <- npar
    class(loglik) <- "logLik"
    return(loglik)
  }
  myidx <- switch(type[1],
                  sim =  seq_along(index(object)$manifest),
                  cond = { endo.idx },
                  cond2 = { endo.idx },
                  exo =  { exo.idx } )

  mom <- moments(object, p, conditional=(type[1]=="cond2"), data=data)
  C <- mom$C
  xi <- mom$xi
  if (type[1]=="exo") {
    C <- C[exo.idx,exo.idx,drop=FALSE]
    xi <- xi[exo.idx,drop=FALSE]
  }
  Debug(list("C=",C),debug)
  k <- nrow(C)

  iC <- Inverse(C,det=TRUE)
  detC <- attributes(iC)$det

  if (!is.null(weight)) {
    weight <- cbind(weight)
    K <- length(exo.idx)+length(endo.idx)
    if (ncol(weight)!=1 & ncol(weight)!=K) {
      w.temp <- weight
      weight <- matrix(1,nrow=nrow(weight),ncol=K)
      weight[,endo.idx] <- w.temp
    }
    if (type=="exo")
      weight <- NULL
  }

  notdatalist <- (!is.list(data) | is.data.frame(data))
  if (missing(n))
    if (!missing(data)) n <- NROW(data)
  if (!missing(n))
  if (notdatalist & (n<2 | indiv | !is.null(weight))) {
    if (n==1)
      data <- rbind(data)
    res <- numeric(n)
    data <- data[,index(object)$manifest,drop=FALSE]
    loglik <- 0;
    for (i in seq_len(n)) {
      ti <- cbind(as.numeric(data[i,myidx]))
      if (meanstructure) {
        ti <- ti-xi
      }
      if (!is.null(weight)) {
        W <- diag(weight[i,],nrow=length(weight[i,]))
        val <- -k/2*log(2*base::pi) -1/2*log(detC) - 1/2*(t(ti)%*%W)%*%iC%*%(ti)
      } else {
        val <- -k/2*log(2*base::pi) -1/2*log(detC) - 1/2*t(ti)%*%iC%*%(ti)
      }
      if (indiv)
        res[i] <- val
      loglik <- loglik + val
    }
    if (indiv)
      return(res)
  } else {
   if (missing(S)) {
      d0 <- procdata.lvm(object,data=data)
      S <- d0$S; mu <- d0$mu; n <- d0$n
    }
    S <- S[myidx,myidx,drop=FALSE]
    mu <- mu[myidx,drop=FALSE]
    T <- S
    if (meanstructure) {
      W <- tcrossprod(mu-xi)
      T <- S+W
    }
    loglik <- -(n*k)/2*log(2*base::pi) -n/2*(log(detC) + tr(T%*%iC))
  }
  return(loglik)
}

###}}}

###{{{ logLik.lvmfit

##' @export
logLik.lvmfit <- function(object,
                          p=coef(object),
                          data=model.frame(object),
                          model=object$estimator,
                          weight=Weight(object),
                          weight2=object$data$weight2,
                          ...) {

  logLikFun <- paste0(model,"_logLik.lvm")
  if (!exists(logLikFun)) {
    model <- "gaussian"
  }
  l <- logLik.lvm(object$model0,p,data,model=model,weight=weight,
                  weight2=weight2,
                  ...)
  return(l)
}

###}}} logLik.lvmfit

###{{{ logLik.lvm.missing

##' @export
logLik.lvm.missing <- function(object,
                               p=pars(object), model=object$estimator,
                               weight=Weight(object$estimate),
                               ...) {
  logLik(object$estimate$model0, p=p, model=model, weight=weight, ...)
}

###}}}

###{{{ logLik.multigroup

##' @export
logLik.multigroup <- function(object,p,data=object$data,weight=NULL,type=c("cond","sim","exo","sat"),...) {
  res <- procrandomslope(object)
  pp <- with(res, modelPar(model,p)$p)

  if (type[1]=="sat") {
    n <- 0
    df <- 0
    loglik <- 0
    for (i in seq_len(object$ngroup)) {
      m <- Model(object)[[i]]
      L <- logLik(m,p=pp[[i]],data=object$data[[i]],type="sat")
      df <- df + attributes(L)$df
      loglik <- loglik + L
      n <- n + object$samplestat[[i]]$n
    }
    attr(loglik, "nall") <- n
    attr(loglik, "nobs") <- n##-df
    attr(loglik, "df") <- df
    class(loglik) <- "logLik"
    return(loglik)
  }

  n <- 0
  loglik <- 0; for (i in seq_len(object$ngroup)) {
    n <- n + object$samplestat[[i]]$n
    val <- logLik(object$lvm[[i]],pp[[i]],data[[i]],weight=weight[[i]],type=type,...)
    loglik <- loglik + val
  }
  attr(loglik, "nall") <- n
  attr(loglik, "nobs") <- n##-length(p)
  attr(loglik, "df") <- length(p)
  class(loglik) <- "logLik"
  return(loglik)
}

###}}} logLik.multigroup

###{{{ logLik.multigroupfit

##' @export
logLik.multigroupfit <- function(object,
                                 p=pars(object), weight=Weight(object), model=object$estimator, ...) {
  logLik(object$model0,p=p,weight=weight,model=model,...)
}
###}}} logLik.multigroup
