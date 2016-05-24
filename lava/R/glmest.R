
glm.estimate.hook <- function(x,estimator,...) {
  yy <- c()
  if (estimator=="glm") {
    for (y in endogenous(x)) {
      fam <- attributes(distribution(x)[[y]])$family
      if (is.null(fam)) fam <- stats::gaussian()
      if (!(tolower(fam$family)%in%
            c("gaussian","gamma","inverse.gaussian"))) {
        yy <- c(yy,y)
      }
    }
    if (length(yy)>0) covariance(x,yy) <- 1
  }
  return(c(list(x=x,estimator=estimator,...)))
}

GLMest <- function(m,data,control=list(),...) {
    v <- vars(m)
    yvar <- endogenous(m)
    res <- c()
    count <- 0
    V <- NULL
    mymsg <- c()
    iids <- c()
    breads <- c()
    for (y in yvar) {
        count <- count+1
        xx <- parents(m,y)
        fam <- attributes(distribution(m)[[y]])$family
        if (is.null(fam)) fam <- stats::gaussian()
        mymsg <- c(mymsg, with(fam, paste0(family,"(",link,")")))
        if (length(xx)==0) xx <- 1
        f <- as.formula(paste0(y,"~",paste(xx,collapse="+")))
        isSurv <- inherits(data[1,y],"Surv")
        if (isSurv) {            
            g <- survival::survreg(f,data=data,dist=fam$family)            
        } else {
            g <- glm(f,family=fam,data=data)
        }
        p <- pars(g)
        ii <- iid(g)
        V0 <- attr(ii,"bread")
        iids <- cbind(iids,ii)
        names(p)[1] <- y
        if (length(p)>1) {
            nn <- paste(y,xx,sep=lava.options()$symbol[1])
            names(p)[seq_along(nn)+1] <- nn
            if (length(p)>length(nn)+1) names(p)[length(p)] <- paste(y,y,sep=lava.options()$symbol[2])
        }
        if (tolower(fam$family)%in%c("gaussian","gamma","inverse.gaussian") && !isSurv) {
            iids <- cbind(iids,0)
            null <- matrix(0); dimnames(null) <- list("scale","scale")
            V0 <- blockdiag(V0,null,pad=0)
        }
        breads <- c(breads,list(V0))
        res <- c(res, list(p));
    }
    coefs <- unlist(res)
    idx <- match(coef(m),names(coefs))
    coefs <- coefs[idx]
    ##V <- Reduce(blockdiag,breads)[idx,idx]
    V <- crossprod(iids[,idx])
    ##V <- crossprod(iids[,idx])    
    mymsg <- noquote(cbind(mymsg))
    colnames(mymsg) <- "Family(Link)"; rownames(mymsg) <- paste(yvar,":")
    list(estimate=coefs,vcov=V,breads=breads,iid=iids[,idx],summary.message=function(...)  {
        mymsg }, dispname="Dispersion:")
}

GLMscore <- function(x,p,data,indiv=TRUE,logLik=FALSE,...) {
    v <- vars(x)
    yvar <- endogenous(x)
    S <- pnames <- c()
    count <- 0
    pos <- 0
    breads <- c()
    L <- 0
    for (y in yvar) {
        count <- count+1
        xx <- parents(x,y)
        pname <- c(y,paste0(y,sep=lava.options()$symbol[1],xx),paste(y,y,sep=lava.options()$symbol[2]))
        pidx <- match(pname,coef(x))
        fam <- attributes(distribution(x)[[y]])$family
        if (is.null(fam)) fam <- stats::gaussian()
        if (length(xx)==0) xx <- 1
        f <- as.formula(paste0(y,"~",paste(xx,collapse="+")))
        isSurv <- inherits(data[1,y],"Surv")
        if (inherits(data[,y],"Surv")) {            
            g <- survival::survreg(f,data=data,dist=fam$family)            
        } else {
            g <- glm(f,family=fam,data=data)
        }
        pdispersion <- NULL
        npar <- length(xx)+2
        p0 <- p[pidx]
        if (!isSurv) L0 <- logL.glm(g,p=p0,indiv=TRUE,...)
        if (tolower(fam$family)%in%c("gaussian","gamma","inverse.gaussian") && !isSurv) {
            p0 <- p0[-length(p0)]
            S0 <- score(g,p=p0,indiv=TRUE,pearson=TRUE,...)
            V0 <- attr(S0,"bread")
            r <- attr(S0,"pearson")
            dispersion <- mean(r^2)
            S0 <- cbind(S0,scale=0)
            null <- matrix(0); dimnames(null) <- list("scale","scale")
            V0 <- blockdiag(V0,null,pad=0)
        } else {
            S0 <- score(g,p=p0,indiv=TRUE,...)
            if (isSurv) L0 <- attr(S0,"logLik")
            V0 <- attr(S0,"bread")
        }
        L <- L+sum(L0)
        breads <- c(breads,list(V0))
        S <- c(S,list(S0))
        pnames <- c(pnames, list(pname));
    }
    coefs <- unlist(pnames)
    idx <- na.omit(match(coef(x),coefs))
    V <- Reduce(blockdiag,breads)[idx,idx]
    S <- Reduce(cbind,S)[,idx,drop=FALSE]
    colnames(S) <- coef(x)
    attributes(S)$bread <- V
    attributes(S)$logLik <- structure(L,nobs=nrow(data),nall=nrow(data),df=length(p),class="logLik")
    if (!indiv) S <- colSums(S)
    return(S)    
}


##' @export
score.lm <- function(x,p=coef(x),data,indiv=FALSE,
                      y,X,offset=NULL,weights=NULL,...) {
  response <- all.vars(formula(x))[1]
  sigma2 <- summary(x)$sigma^2
  if (missing(data)) {
      X <- model.matrix(x)
      y <- model.frame(x)[,1]
  } else {
      X <- model.matrix(formula(x),data=data)
      y <- model.frame(formula(x),data=data)[,1]
  }
  n <- nrow(X)
  if(any(is.na(p))) warning("Over-parameterized model")
  Xbeta <- X%*%p
  if (is.null(offset)) offset <- x$offset
  if (!is.null(offset)) Xbeta <- Xbeta+offset
  r <- y-Xbeta
  if (is.null(weights)) weights <- x$weights
  if (!is.null(weights)) r <- r*weights
  A <- as.vector(r)/sigma2
  S <- apply(X,2,function(x) x*A)
  if (!indiv) return(colSums(S))
  attributes(S)$bread <- vcov(x)
  return(S)
}

##' @export
score.glm <- function(x,p=coef(x),data,indiv=FALSE,pearson=FALSE,
                      y,X,link,dispersion,offset=NULL,weights=NULL,...) {

    response <- all.vars(formula(x))[1]
    if (inherits(x,"glm")) {
        link <- family(x)
        if (missing(data)) {
            X <- model.matrix(x)
            y <- model.frame(x)[,1]
        } else {
            X <- model.matrix(formula(x),data=data)
            y <- model.frame(formula(x),data=data)[,1]
        }
        offset <- x$offset
    } else {
        if (missing(link)) stop("Family needed")
        if (missing(data)) stop("data needed")
        X <- model.matrix(formula(x),data=data)
        y <- model.frame(formula(x),data=data)[,1]
    }
    if (is.character(y) || is.factor(y)) {
        y <- as.numeric(as.factor(y))-1
    }
    n <- nrow(X)
    g <- link$linkfun
    ginv <- link$linkinv
    dginv <- link$mu.eta ## D[linkinv]
    ##dg <- function(x) 1/dginv(g(x)) ## Dh^-1 = 1/(h'(h^-1(x)))
    canonf <- do.call(link$family,list())
    caninvlink <- canonf$linkinv
    canlink <- canonf$linkfun
    Dcaninvlink <- canonf$mu.eta
    Dcanlink <- function(x) 1/Dcaninvlink(canlink(x))
    ##gmu <- function(x) g(caninvlink(x))
    ##invgmu <- function(z) canlink(ginv(z))
    h <- function(z) Dcanlink(ginv(z))*dginv(z)
    if(any(is.na(p))) stop("Over-parameterized model")
    Xbeta <- X%*%p
    if (!is.null(offset)) Xbeta <- Xbeta+offset
    if (missing(data) && !is.null(x$offset) && is.null(offset) ) Xbeta <- Xbeta+x$offset
    pi <- ginv(Xbeta)
    ##res <- as.vector(y/pi*dginv(Xbeta)-(1-y)/(1-pi)*dginv(Xbeta))*X
    ##return(res)
    r <- y-pi
    if (!is.null(x$prior.weights) || !is.null(weights)) {
        if (is.null(weights)) weights <- x$prior.weights
    } else {
        weights <- !is.na(r)
    }
    r <- r*weights
    a.phi <- 1
    rpearson <- as.vector(r)/link$variance(pi)^.5
    if (length(p)>length(coef(x))) {
        a.phi <- p[length(coef(x))+1]
    } else if (tolower(family(x)$family)%in%c("gaussian","gamma","inverse.gaussian")) {
        ##a.phi <- summary(x)$dispersion*g0$df.residual/sum(weights)
        a.phi <- sum(rpearson^2)*x$df.residual/x$df.residual^2
    }
    A <- as.vector(h(Xbeta)*r)/a.phi
    S <- apply(X,2,function(x) x*A)
    if (!indiv) return(colSums(S))
    if (pearson) attr(S,"pearson") <- rpearson
    attributes(S)$bread <- vcov(x)
    if (x$family$family=="quasi" && x$family$link=="identity" && x$family$varfun=="constant")
        attributes(S)$bread <- -Inverse(information.glm(x))
    return(S)
}

##' @export
pars.glm <- function(x,...) {
  if (tolower(family(x)$family)%in%c("gaussian","gamma","inverse.gaussian")) {
    res <- c(coef(x),summary(x)$dispersion)
    names(res)[length(res)] <- "Dispersion"
    return(res)
  }
  return(coef(x))
}

logL.glm <- function(x,p=pars.glm(x),data,indiv=FALSE,...) {
    if (!missing(data)) {
        x <- update(x,data=data,...)
    }
    f <- family(x)
    ginv <- f$linkinv
    X <- model.matrix(x)
    n <- nrow(X)
    disp <- 1; p0 <- p
    if (tolower(family(x)$family)%in%c("gaussian","gamma","inverse.gaussian")) {
        if (length(p)==ncol(X)) {
            disp <- suppressWarnings((summary(x)$dispersion))
        } else {
            disp <- tail(p,1)
            p0 <- p[-length(p)]
        }
    }
    if(any(is.na(p))) stop("Over-parametrized model")
    Xbeta <- X%*%p0
    if (!is.null(x$offset)) Xbeta <- Xbeta+x$offset
    y <- model.frame(x)[,1]
    mu <- ginv(Xbeta)
    w <- x$prior.weights
    dev <-  f$dev.resids(y,mu,w)
    if (indiv) {

    }
    loglik <- length(p)-(f$aic(y,n,mu,w,sum(dev))/2+x$rank)
    structure(loglik,nobs=n,df=length(p),class="logLik")
}

##' @export
iid.glm <- function(x,...) {
    ## if (x$family$family=="quasi" && x$family$link=="identity" && x$family$varfun=="constant") {
    ##     return(iid.default(x,information.glm,...))
    ## }
    iid.default(x,...)
}

hessian.glm <- function(x,p=coef(x),...) {
  numDeriv::jacobian(function(theta) score.glm(x,p=theta,indiv=FALSE,...),p)
}

##' @export
information.glm <- function(x,...) hessian.glm(x,...)

robustvar <- function(x,id=NULL,...) {
  U <- score(x,indiv=TRUE)
  II <- unique(id)
  K <- length(II)
  J <- 0
  if (is.null(id)) {
    J <- crossprod(U)
  } else {
    for (ii in II) {
      J <- J+tcrossprod(colSums(U[which(id==ii),,drop=FALSE]))
    }
    J <- K/(K-1)*J
  }
  iI <- vcov(x)
  V <- iI%*%J%*%iI
  return(V)
}

glm_logLik.lvm <- function(object,...) {
    attr(GLMscore(object,...),"logLik")
}

glm_method.lvm <- NULL
glm_objective.lvm <- function(x,p,data,...) {
  GLMest(x,data,...)
}
glm_gradient.lvm <- function(x,p,data,...) {
  -GLMscore(x,p,data,...)
}

glm_variance.lvm <- function(x,p,data,opt,...) {
  opt$vcov
}
