## Fit gssanova0 model
gssanova0 <- function(formula,family,type=NULL,data=list(),weights,
                      subset,offset,na.action=na.omit,partial=NULL,
                      method=NULL,varht=1,nu=NULL,prec=1e-7,maxiter=30)
{
    ## Obtain model frame and model terms
    mf <- match.call()
    mf$family <- mf$type <- mf$partial <- NULL
    mf$method <- mf$varht <- mf$nu <- NULL
    mf$prec <- mf$maxiter <- NULL
    mf[[1]] <- as.name("model.frame")
    mf <- eval(mf,parent.frame())
    ## Generate terms
    term <- mkterm(mf,type)
    ## Specify default method
    if (is.null(method)) {
        method <- switch(family,
                         binomial="u",
                         nbinomial="u",
                         poisson="u",
                         inverse.gaussian="v",
                         Gamma="v",
                         weibull="u",
                         lognorm="u",
                         loglogis="u")
    }
    ## Generate s, q, and y
    nobs <- dim(mf)[1]
    s <- q <- NULL
    nq <- 0
    for (label in term$labels) {
        if (label=="1") {
            s <- cbind(s,rep(1,len=nobs))
            next
        }
        x <- mf[,term[[label]]$vlist]
        nphi <- term[[label]]$nphi
        nrk <- term[[label]]$nrk
        if (nphi) {
            phi <- term[[label]]$phi
            for (i in 1:nphi)
                s <- cbind(s,phi$fun(x,nu=i,env=phi$env))
        }
        if (nrk) {
            rk <- term[[label]]$rk
            for (i in 1:nrk) {
                nq <- nq+1
                q <- array(c(q,rk$fun(x,x,nu=i,env=rk$env,out=TRUE)),c(nobs,nobs,nq))
            }
        }
    }
    ## Add the partial term
    if (!is.null(partial)) {
        mf.p <- model.frame(partial,data)
        for (lab in colnames(mf.p)) mf[,lab] <- mf.p[,lab]
        mt.p <- attr(mf.p,"terms")
        lab.p <- labels(mt.p)
        matx.p <- model.matrix(mt.p,data)[,-1,drop=FALSE]
        if (dim(matx.p)[1]!=dim(mf)[1])
            stop("gss error in ssanova: partial data are of wrong size")
        matx.p <- scale(matx.p)
        center.p <- attr(matx.p,"scaled:center")
        scale.p <- attr(matx.p,"scaled:scale")
        s <- cbind(s,matx.p)
        part <- list(mt=mt.p,center=center.p,scale=scale.p)
    }
    else part <- lab.p <- NULL
    if (qr(s)$rank<dim(s)[2])
        stop("gss error in gssanova0: unpenalized terms are linearly dependent")
    y <- model.response(mf,"numeric")
    wt <- model.weights(mf)
    offset <- model.offset(mf)
    if (!is.null(offset)) {
        term$labels <- c(term$labels,"offset")
        term$offset <- list(nphi=0,nrk=0)
    }
    if (!nq) stop("gss error in gssanova0: use glm for models with only unpenalized terms")
    ## Fit the model
    if (nq==1) {
        q <- q[,,1]
        z <- sspregpoi(family,s,q,y,wt,offset,method,varht,nu,prec,maxiter)
    }
    else z <- mspregpoi(family,s,q,y,wt,offset,method,varht,nu,prec,maxiter)
    ## Brief description of model terms
    desc <- NULL
    for (label in term$labels)
        desc <- rbind(desc,as.numeric(c(term[[label]][c("nphi","nrk")])))
    if (!is.null(partial)) {
        desc <- rbind(desc,matrix(c(1,0),length(lab.p),2,byrow=TRUE))
    }
    desc <- rbind(desc,apply(desc,2,sum))
    if (is.null(partial)) rownames(desc) <- c(term$labels,"total")
    else rownames(desc) <- c(term$labels,lab.p,"total")
    colnames(desc) <- c("Unpenalized","Penalized")
    ## Return the results
    obj <- c(list(call=match.call(),family=family,mf=mf,terms=term,desc=desc,
                  partial=part,lab.p=lab.p),z)
    class(obj) <- c("gssanova0","ssanova0","gssanova")
    obj
}

## Fit Single Smoothing Parameter REGression by Performance-Oriented Iteration
sspregpoi <- function(family,s,q,y,wt,offset,method="u",
                      varht=1,nu,prec=1e-7,maxiter=30)
{
    ## Check inputs
    if (is.vector(s)) s <- as.matrix(s)
    if (!(is.matrix(s)&is.matrix(q)&is.character(method))) {
        stop("gss error in sspregpoi: inputs are of wrong types")
    }
    nobs <- dim(s)[1]
    nnull <- dim(s)[2]
    if (!((dim(s)[1]==nobs)&(dim(q)[1]==nobs)&(dim(q)[2]==nobs)
          &(nobs>=nnull)&(nnull>0))) {
        stop("gss error in sspregpoi: inputs have wrong dimensions")
    }
    ## Set method for smoothing parameter selection
    code <- (1:3)[c("v","m","u")==method]
    if (!length(code)) {
        stop("gss error: unsupported method for smoothing parameter selection")
    }
    eta <- rep(0,nobs)
    nla0 <- log10(mean(abs(diag(q))))
    limnla <- nla0+c(-.5,.5)
    iter <- 0
    if (family=="nbinomial") nu <- NULL
    else nu <- list(nu,is.null(nu))
    repeat {
        iter <- iter+1
        dat <- switch(family,
                      binomial=mkdata.binomial(y,eta,wt,offset),
                      nbinomial=mkdata.nbinomial(y,eta,wt,offset,nu),
                      poisson=mkdata.poisson(y,eta,wt,offset),
                      inverse.gaussian=mkdata.inverse.gaussian(y,eta,wt,offset),
                      Gamma=mkdata.Gamma(y,eta,wt,offset),
                      weibull=mkdata.weibull(y,eta,wt,offset,nu),
                      lognorm=mkdata.lognorm(y,eta,wt,offset,nu),
                      loglogis=mkdata.loglogis(y,eta,wt,offset,nu))
        nu <- dat$nu
        w <- as.vector(sqrt(dat$wt))
        ywk <- w*dat$ywk
        swk <- w*s
        qwk <- w*t(w*q)
        ## Call RKPACK driver DSIDR
        z <- .Fortran("dsidr0",
                      as.integer(code),
                      swk=as.double(swk), as.integer(nobs),
                      as.integer(nobs), as.integer(nnull),
                      as.double(ywk),
                      qwk=as.double(qwk), as.integer(nobs),
                      as.double(0), as.integer(-1), as.double(limnla),
                      nlambda=double(1), score=double(1), varht=as.double(varht),
                      c=double(nobs), d=double(nnull),
                      qraux=double(nnull), jpvt=integer(nnull),
                      double(3*nobs),
                      info=integer(1),PACKAGE="gss")
        ## Check info for error
        if (info<-z$info) {               
            if (info>0)
                stop("gss error in sspregpoi: matrix s is rank deficient")
            if (info==-2)
                stop("gss error in sspregpoi: matrix q is indefinite")
            if (info==-1)
                stop("gss error in sspregpoi: input data have wrong dimensions")
            if (info==-3)
                stop("gss error in sspregpoi: unknown method for smoothing parameter selection.")
        }
        eta.new <- (ywk-10^z$nlambda*z$c)/w
        if (!is.null(offset)) eta.new <- eta.new + offset
        disc <- sum(dat$wt*((eta-eta.new)/(1+abs(eta)))^2)/sum(dat$wt)
        limnla <- pmax(z$nlambda+c(-.5,.5),nla0-5)
        if (disc<prec) break
        if (iter>=maxiter) {
            warning("gss warning in gssanova0: performance-oriented iteration fails to converge")
            break
        }
        eta <- eta.new
    }
    ## Return the fit
    if (is.list(nu)) nu <- nu[[1]]
    c(list(method=method,theta=0,w=as.vector(dat$wt),
           eta=as.vector(eta),iter=iter,nu=nu),
      z[c("c","d","nlambda","score","varht","swk","qraux","jpvt","qwk")])
}

## Fit Multiple Smoothing Parameter REGression by Performance-Oriented Iteration
mspregpoi <- function(family,s,q,y,wt,offset,method="u",
                      varht=1,nu,prec=1e-7,maxiter=30)
{
    ## Check inputs
    if (is.vector(s)) s <- as.matrix(s)
    if (!(is.matrix(s)&is.array(q)&(length(dim(q))==3)
          &is.character(method))) {
        stop("gss error in mspregpoi: inputs are of wrong types")
    }
    nobs <- dim(s)[1]
    nnull <- dim(s)[2]
    nq <- dim(q)[3]
    if (!((dim(s)[1]==nobs)&(dim(q)[1]==nobs)&(dim(q)[2]==nobs)
          &(nobs>=nnull)&(nnull>0))) {
        stop("gss error in sspregpoi: inputs have wrong dimensions")
    }
    ## Set method for smoothing parameter selection
    code <- (1:3)[c("v","m","u")==method]
    if (!length(code)) {
        stop("gss error: unsupported method for smoothing parameter selection")
    }
    eta <- rep(0,nobs)
    init <- 0
    theta <- rep(0,nq)
    iter <- 0
    if (family=="nbinomial") nu <- NULL
    else nu <- list(nu,is.null(nu))
    qwk <- array(0,c(nobs,nobs,nq))
    repeat {
        iter <- iter+1
        dat <- switch(family,
                      binomial=mkdata.binomial(y,eta,wt,offset),
                      nbinomial=mkdata.nbinomial(y,eta,wt,offset,nu),
                      poisson=mkdata.poisson(y,eta,wt,offset),
                      inverse.gaussian=mkdata.inverse.gaussian(y,eta,wt,offset),
                      Gamma=mkdata.Gamma(y,eta,wt,offset),
                      weibull=mkdata.weibull(y,eta,wt,offset,nu),
                      lognorm=mkdata.lognorm(y,eta,wt,offset,nu),
                      loglogis=mkdata.loglogis(y,eta,wt,offset,nu))
        nu <- dat$nu
        w <- as.vector(sqrt(dat$wt))
        ywk <- w*dat$ywk
        swk <- w*s
        for (i in 1:nq) qwk[,,i] <- w*t(w*q[,,i])
        ## Call RKPACK driver DMUDR
        z <- .Fortran("dmudr0",
                      as.integer(code),
                      as.double(swk),   # s
                      as.integer(nobs), as.integer(nobs), as.integer(nnull),
                      as.double(qwk),   # q
                      as.integer(nobs), as.integer(nobs), as.integer(nq),
                      as.double(ywk),   # y
                      as.double(0), as.integer(init),
                      as.double(prec), as.integer(maxiter),
                      theta=as.double(theta), nlambda=double(1),
                      score=double(1), varht=as.double(varht),
                      c=double(nobs), d=double(nnull),
                      double(nobs*nobs*(nq+2)),
                      info=integer(1),PACKAGE="gss")[c("theta","nlambda","c","info")]
        ## Check info for error
        if (info<-z$info) {               
            if (info>0)
                stop("gss error in mspreg: matrix s is rank deficient")
            if (info==-2)
                stop("gss error in mspreg: matrix q is indefinite")
            if (info==-1)
                stop("gss error in mspreg: input data have wrong dimensions")
            if (info==-3)
                stop("gss error in mspreg: unknown method for smoothing parameter selection.")
            if (info==-4)
                stop("gss error in mspreg: iteration fails to converge, try to increase maxiter")
            if (info==-5)
                stop("gss error in mspreg: iteration fails to find a reasonable descent direction")
        }
        eta.new <- (ywk-10^z$nlambda*z$c)/w
        if (!is.null(offset)) eta.new <- eta.new + offset
        disc <- sum(dat$wt*((eta-eta.new)/(1+abs(eta)))^2)/sum(dat$wt)
        if (disc<prec) break
        if (iter>=maxiter) {
            warning("gss warning in gssanova0: performance-oriented iteration fails to converge")
            break
        }
        init <- 1
        theta <- z$theta
        eta <- eta.new
    }
    qqwk <- 10^z$theta[1]*qwk[,,1]
    for (i in 2:nq) qqwk <- qqwk + 10^z$theta[i]*qwk[,,i]
    ## Call RKPACK driver DSIDR
    z <- .Fortran("dsidr0",
                  as.integer(code),
                  swk=as.double(swk), as.integer(nobs),
                  as.integer(nobs), as.integer(nnull),
                  as.double(ywk),
                  qwk=as.double(qqwk), as.integer(nobs),
                  as.double(0), as.integer(0), double(2),
                  nlambda=double(1), score=double(1), varht=as.double(varht),
                  c=double(nobs), d=double(nnull),
                  qraux=double(nnull), jpvt=integer(nnull),
                  double(3*nobs),
                  info=integer(1),PACKAGE="gss")
    ## Check info for error
    if (info<-z$info) {               
        if (info>0)
            stop("gss error in sspregpoi: matrix s is rank deficient")
        if (info==-2)
            stop("gss error in sspregpoi: matrix q is indefinite")
        if (info==-1)
            stop("gss error in sspregpoi: input data have wrong dimensions")
        if (info==-3)
            stop("gss error in sspregpoi: unknown method for smoothing parameter selection.")
    }
    ## Return the fit
    if (is.list(nu)) nu <- nu[[1]]
    c(list(method=method,theta=theta,w=as.vector(dat$wt),
           eta=as.vector(eta),iter=iter,nu=nu),
      z[c("c","d","nlambda","score","varht","swk","qraux","jpvt","qwk")])
}
