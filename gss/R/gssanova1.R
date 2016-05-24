## Fit gssanova model
gssanova1 <- function(formula,family,type=NULL,data=list(),weights,
                      subset,offset,na.action=na.omit,partial=NULL,
                      method=NULL,varht=1,alpha=1.4,nu=NULL,
                      id.basis=NULL,nbasis=NULL,seed=NULL,random=NULL,
                      skip.iter=FALSE)
{
    if (!(family%in%c("binomial","poisson","Gamma","nbinomial","inverse.gaussian",
                      "weibull","lognorm","loglogis")))
        stop("gss error in gssanova1: family not implemented")
    ## Obtain model frame and model terms
    mf <- match.call()
    mf$family <- mf$type <- mf$partial <- NULL
    mf$method <- mf$varht <- mf$nu <- NULL
    mf$alpha <- mf$id.basis <- mf$nbasis <- mf$seed <- NULL
    mf$random <- mf$skip.iter <- NULL
    mf[[1]] <- as.name("model.frame")
    mf <- eval(mf,parent.frame())
    wt <- model.weights(mf)
    ## Generate sub-basis
    nobs <- dim(mf)[1]
    if (is.null(id.basis)) {
        if (is.null(nbasis))  nbasis <- max(30,ceiling(10*nobs^(2/9)))
        if (nbasis>=nobs)  nbasis <- nobs
        if (!is.null(seed))  set.seed(seed)
        id.basis <- sample(nobs,nbasis,prob=wt)
    }
    else {
        if (max(id.basis)>nobs|min(id.basis)<1)
            stop("gss error in gssanova: id.basis out of range")
        nbasis <- length(id.basis)
    }
    ## Generate terms
    term <- mkterm(mf,type)
    ## Generate random
    if (!is.null(random)) {
        if (class(random)=="formula") random <- mkran(random,data)
    }
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
    ## Generate s, r, and y
    s <- r <- NULL
    nq <- 0
    for (label in term$labels) {
        if (label=="1") {
            s <- cbind(s,rep(1,len=nobs))
            next
        }
        x <- mf[,term[[label]]$vlist]
        x.basis <- mf[id.basis,term[[label]]$vlist]
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
                r <- array(c(r,rk$fun(x,x.basis,nu=i,env=rk$env,out=TRUE)),c(nobs,nbasis,nq))
            }
        }
    }
    if (is.null(r))
        stop("gss error in gssanova1: use glm for models with only unpenalized terms")
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
        stop("gss error in gssanova: unpenalized terms are linearly dependent")
    ## Prepare the data
    y <- model.response(mf,"numeric")
    offset <- model.offset(mf)
    if (!is.null(offset)) {
        term$labels <- c(term$labels,"offset")
        term$offset <- list(nphi=0,nrk=0)
    }
    nu.wk <- list(NULL,FALSE)
    if ((family=="nbinomial")&is.vector(y)) nu.wk <- list(NULL,TRUE)
    if (family%in%c("weibull","lognorm","loglogis")) {
        if (is.null(nu)) nu.wk <- list(nu,TRUE)
        else nu.wk <- list(nu,FALSE)
    }
    ## Fit the model
    z <- ngreg1(family,s,r,id.basis,y,wt,offset,method,varht,alpha,nu.wk,
                random,skip.iter)
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
                  alpha=alpha,id.basis=id.basis,partial=part,lab.p=lab.p,
                  random=random,skip.iter=skip.iter),z)
    class(obj) <- c("gssanova","ssanova")
    obj
}

## Fit Single Smoothing Parameter REGression by Performance-Oriented Iteration
ngreg1 <- function(family,s,r,id.basis,y,wt,offset,method,varht,alpha,nu,
                   random,skip.iter)
{
    nobs <- dim(s)[1]
    eta <- rep(0,nobs)
    iter <- 0
    if (nu[[2]]&is.null(nu[[1]])) {
        eta <- rep(0,nobs)
        wk <- switch(family,
                      nbinomial=mkdata.nbinomial(y,eta,wt,offset,NULL),
                      weibull=mkdata.weibull(y,eta,wt,offset,nu),
                      lognorm=mkdata.lognorm(y,eta,wt,offset,nu),
                      loglogis=mkdata.loglogis(y,eta,wt,offset,nu))
        nu <- wk$nu
    }
    nq <- dim(r)[3]
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
        if (nq==1) {
            z <- sspreg1(s,r[,,1],r[id.basis,,1],dat$ywk,w,method,alpha,varht,random)
        }
        else z <- mspreg1(s,r,id.basis,dat$ywk,w,method,alpha,varht,random,skip.iter)
        r.wk <- 0
        for (i in 1:nq) r.wk <- r.wk + 10^z$theta[i]*r[,,i]
        if (!is.null(random)) r.wk <- cbind(r.wk,random$z)
        eta.new <- as.vector(s%*%z$d+r.wk%*%c(z$c,z$b))
        if (!is.null(offset)) eta.new <- eta.new + offset
        disc <- sum(dat$wt*((eta-eta.new)/(1+abs(eta)))^2)/sum(dat$wt)
        eta <- eta.new
        if (disc<1e-7) break
        if (iter>=30) {
            warning("gss warning in gssanova1: performance-oriented iteration fails to converge")
            break
        }
    }
    ## Return the fit
    if (is.list(nu)) nu <- nu[[1]]
    c(z,list(nu=nu,eta=eta,w=dat$wt))
}
