## Fit gssanova model
gssanova <- function(formula,family,type=NULL,data=list(),weights,
                     subset,offset,na.action=na.omit,partial=NULL,
                     alpha=NULL,nu=NULL,
                     id.basis=NULL,nbasis=NULL,seed=NULL,random=NULL,
                     skip.iter=FALSE)
{
    if (!(family%in%c("binomial","poisson","Gamma","inverse.gaussian","nbinomial",
                      "weibull","lognorm","loglogis")))
        stop("gss error in gssanova: family not implemented")
    if (is.null(alpha)) {
        alpha <- 1.4
        if (family%in%c("binomial","nbinomial","inverse.gaussian")) alpha <- 1
    }
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
        stop("gss error in gssanova: use glm for models with only unpenalized terms")
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
    if (nq==1) {
        r <- r[,,1]
        z <- sspngreg(family,s,r,r[id.basis,],y,wt,offset,alpha,nu.wk,random)
    }
    else z <- mspngreg(family,s,r,id.basis,y,wt,offset,alpha,nu.wk,random,skip.iter)
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

## Fit Single Smoothing Parameter Non-Gaussian REGression
sspngreg <- function(family,s,r,q,y,wt,offset,alpha,nu,random)
{
    nobs <- nrow(r)
    nxi <- ncol(r)
    if (!is.null(s)) {
        if (is.vector(s)) nnull <- 1
        else nnull <- ncol(s)
    }
    else nnull <- 0
    if (!is.null(random)) nz <- ncol(as.matrix(random$z))
    else nz <- 0
    nxiz <- nxi + nz
    nn <- nxiz + nnull
    ## cv function
    cv <- function(lambda) {
        if (nu[[2]]) {
            la.wk <- lambda[-2]
            nu.wk <- list(exp(lambda[2]),FALSE)
        }
        else {
            la.wk <- lambda
            nu.wk <- nu
        }
        if (is.null(random)) q.wk <- 10^(la.wk+theta)*q
        else {
            q.wk <- matrix(0,nxiz,nxiz)
            q.wk[1:nxi,1:nxi] <- 10^(la.wk[1]+theta)*q
            q.wk[(nxi+1):nxiz,(nxi+1):nxiz] <-
                10^(2*ran.scal)*random$sigma$fun(la.wk[-1],random$sigma$env)
        }
        alpha.wk <- max(0,log.la0-la.wk[1]-5)*(3-alpha) + alpha
        alpha.wk <- min(alpha.wk,3)
        z <- ngreg(dc,family,cbind(s,10^theta*r),q.wk,y,wt,offset,nu.wk,alpha.wk)
        assign("dc",z$dc,inherits=TRUE)
        assign("fit",z[c(1:3,5:10)],inherits=TRUE)
        z$score
    }
    cv.wk <- function(lambda) cv.scale*cv(lambda)+cv.shift
    ## initialization
    tmp <- sum(r^2)
    if (is.null(s)) theta <- 0
    else theta <- log10(sum(s^2)/nnull/tmp*nxi) / 2
    log.la0 <- log10(tmp/sum(diag(q))) + theta
    if (!is.null(random)) {
        ran.scal <- theta - log10(sum(random$z^2)/nz/tmp*nxi) / 2
        r <- cbind(r,10^(ran.scal-theta)*random$z)
    }
    else ran.scal <- NULL
    if (nu[[2]]&is.null(nu[[1]])) {
        eta <- rep(0,nobs)
        wk <- switch(family,
                      nbinomial=mkdata.nbinomial(y,eta,wt,offset,NULL),
                      weibull=mkdata.weibull(y,eta,wt,offset,nu),
                      lognorm=mkdata.lognorm(y,eta,wt,offset,nu),
                      loglogis=mkdata.loglogis(y,eta,wt,offset,nu))
        nu[[1]] <- wk$nu[[1]]
    }
    ## lambda search
    dc <- rep(0,nn)
    fit <- NULL
    la <- log.la0
    if (nu[[2]]) la <- c(la, log(nu[[1]]))
    if (!is.null(random)) la <- c(la,random$init)
    if (length(la)-1) {
        counter <- 0
        ## scale and shift cv
        tmp <- abs(cv(la))
        cv.scale <- 1
        cv.shift <- 0
        if (tmp<1&tmp>10^(-4)) {
            cv.scale <- 10/tmp
            cv.shift <- 0
        }
        if (tmp<10^(-4)) {
            cv.scale <- 10^2
            cv.shift <- 10
        }
        repeat {
            zz <- nlm(cv.wk,la,stepmax=1,ndigit=7)
            if (zz$code<=3) break
            la <- zz$est
            counter <- counter + 1
            if (counter>=5) {
                warning("gss warning in ssanova: iteration for model selection fails to converge")
                break
            }
        }
    }
    else {
        mn0 <- log.la0-6
        mx0 <- log.la0+6
        repeat {
            mn <- max(la-1,mn0)
            mx <- min(la+1,mx0)
            zz <- nlm0(cv,c(mn,mx))
            if ((min(zz$est-mn,mx-zz$est)>=1e-1)||
                (min(zz$est-mn0,mx0-zz$est)<1e-1)) break
            else la <- zz$est
        }
    }
    ## return
    jk <- cv(zz$est)
    if (nu[[2]]) {
        nu.wk <- exp(zz$est[2])
        zz$est <- zz$est[-2]
    }
    else nu.wk <- NULL
    if (is.null(random)) q.wk <- 10^theta*q
    else {
        q.wk <- matrix(0,nxiz,nxiz)
        q.wk[1:nxi,1:nxi] <- 10^theta*q
        q.wk[(nxi+1):nxiz,(nxi+1):nxiz] <-
            10^(2*ran.scal-zz$est[1])*random$sigma$fun(zz$est[-1],random$sigma$env)
    }
    se.aux <- regaux(sqrt(fit$w)*s,10^theta*sqrt(fit$w)*r,q.wk,zz$est[1],fit)
    c <- fit$dc[nnull+(1:nxi)]
    if (nnull) d <- fit$dc[1:nnull]
    else d <- NULL
    if (nz) b <- 10^(ran.scal)*fit$dc[nnull+nxi+(1:nz)]
    else b <- NULL
    c(list(theta=theta,ran.scal=ran.scal,c=c,d=d,b=b,nlambda=zz$est[1],
           zeta=zz$est[-1],nu=nu.wk),fit[-1],list(se.aux=se.aux))
}

## Fit Multiple Smoothing Parameter Non-Gaussian REGression
mspngreg <- function(family,s,r,id.basis,y,wt,offset,alpha,nu,random,skip.iter)
{
    nobs <- nrow(r)
    nxi <- ncol(r)
    if (!is.null(s)) {
        if (is.vector(s)) nnull <- 1
        else nnull <- ncol(s)
    }
    else nnull <- 0
    if (!is.null(random)) nz <-ncol(as.matrix(random$z))
    else nz <- 0
    nxiz <- nxi + nz
    nn <- nxiz + nnull
    nq <- dim(r)[3]
    ## cv function
    cv <- function(theta) {
        if (nu[[2]]) {
            the.wk <- theta[-(nq+1)]
            nu.wk <- list(exp(theta[nq+1]),FALSE)
        }
        else {
            the.wk <- theta
            nu.wk <- nu
        }
        ind.wk <- theta[1:nq]!=theta.old
        if (sum(ind.wk)==nq) {
            r.wk0 <- 0
            for (i in 1:nq) {
                r.wk0 <- r.wk0 + 10^theta[i]*r[,,i]
            }
            assign("r.wk",r.wk0+0,inherits=TRUE)
            assign("theta.old",theta[1:nq]+0,inherits=TRUE)
        }
        else {
            r.wk0 <- r.wk
            for (i in (1:nq)[ind.wk]) {
                theta.wk <- (10^(theta[i]-theta.old[i])-1)*10^theta.old[i]
                r.wk0 <- r.wk0 + theta.wk*r[,,i]
            }
        }
        qq.wk <- r.wk0[id.basis,]
        if (is.null(random)) q.wk <- 10^nlambda*qq.wk
        else {
            r.wk0 <- cbind(r.wk0,10^(ran.scal)*random$z)
            q.wk <- matrix(0,nxiz,nxiz)
            q.wk[1:nxi,1:nxi] <- 10^nlambda*qq.wk
            q.wk[(nxi+1):nxiz,(nxi+1):nxiz] <-
                10^(2*ran.scal)*random$sigma$fun(the.wk[-(1:nq)],random$sigma$env)
        }
        alpha.wk <- max(0,the.wk[1:nq]-log.th0-5)*(3-alpha) + alpha
        alpha.wk <- min(alpha.wk,3)
        z <- ngreg(dc,family,cbind(s,r.wk0),q.wk,y,wt,offset,nu.wk,alpha.wk)
        assign("dc",z$dc,inherits=TRUE)
        assign("fit",z[c(1:3,5:10)],inherits=TRUE)
        z$score
    }
    cv.wk <- function(theta) cv.scale*cv(theta)+cv.shift
    ## initialization
    theta <- -log10(apply(r[id.basis,,],3,function(x)sum(diag(x))))
    r.wk <- 0
    for (i in 1:nq) {
        r.wk <- r.wk + 10^theta[i]*r[,,i]
    }
    ## theta adjustment
    z <- sspngreg(family,s,r.wk,r.wk[id.basis,],y,wt,offset,alpha,nu,random)
    if (nu[[2]]) nu[[1]] <- z$nu
    theta <- theta + z$theta
    r.wk <- 0
    for (i in 1:nq) {
        theta[i] <- 2*theta[i] + log10(t(z$c)%*%r[id.basis,,i]%*%z$c)
        r.wk <- r.wk + 10^theta[i]*r[,,i]
    }
    log.la0 <- log10(sum(r.wk^2)/sum(diag(r.wk[id.basis,])))
    log.th0 <- theta-log.la0
    ## lambda search
    z <- sspngreg(family,s,r.wk,r.wk[id.basis,],y,wt,offset,alpha,nu,random)
    if (nu[[2]]) nu[[1]] <- z$nu
    nlambda <- z$nlambda
    log.th0 <- log.th0 + z$lambda
    theta <- theta + z$theta
    if (!is.null(random)) ran.scal <- z$ran.scal
    ## early return
    if (skip.iter) {
        z$theta <- theta
        return(z)
    }
    ## theta search
    dc <- rep(0,nn)
    fit <- NULL
    theta.old <- theta
    if (nu[[2]]) theta <- c(theta, log(nu[[1]]))
    if (!is.null(random)) theta <- c(theta,z$zeta)
    counter <- 0
    r.wk <- 0
    for (i in 1:nq) {
        r.wk <- r.wk + 10^theta[i]*r[,,i]
    }
    tmp <- abs(cv(theta))
    cv.scale <- 1
    cv.shift <- 0
    if (tmp<1&tmp>10^(-4)) {
        cv.scale <- 10/tmp
        cv.shift <- 0
    }
    if (tmp<10^(-4)) {
        cv.scale <- 10^2
        cv.shift <- 10
    }
    repeat {
        zz <- nlm(cv.wk,theta,stepmax=1,ndigit=7)
        if (zz$code<=3)  break
        theta <- zz$est        
        counter <- counter + 1
        if (counter>=5) {
            warning("gss warning in gssanova: iteration for model selection fails to converge")
            break
        }
    }
    ## return
    jk <- cv(zz$est)
    if (nu[[2]]) {
        nu.wk <- exp(zz$est[nq+1])
        zz$est <- zz$est[-(nq+1)]
    }
    else nu.wk <- NULL
    r.wk <- 0
    for (i in 1:nq) {
        r.wk <- r.wk + 10^zz$est[i]*r[,,i]
    }
    qq.wk <- r.wk[id.basis,]
    if (is.null(random)) q.wk <- qq.wk
    else {
        r.wk <- cbind(r.wk,10^(ran.scal)*random$z)
        q.wk <- matrix(0,nxiz,nxiz)
        q.wk[1:nxi,1:nxi] <- qq.wk
        q.wk[(nxi+1):nxiz,(nxi+1):nxiz] <-
            10^(2*ran.scal-nlambda)*random$sigma$fun(zz$est[-(1:nq)],random$sigma$env)
    }
    se.aux <- regaux(sqrt(fit$w)*s,sqrt(fit$w)*r.wk,q.wk,nlambda,fit)
    c <- fit$dc[nnull+(1:nxi)]
    if (nnull) d <- fit$dc[1:nnull]
    else d <- NULL
    if (nz) b <- 10^(ran.scal)*fit$dc[nnull+nxi+(1:nz)]
    else b <- NULL
    c(list(theta=zz$est[1:nq],c=c,d=d,b=b,nlambda=nlambda,zeta=zz$est[-(1:nq)],nu=nu.wk),
      fit[-1],list(se.aux=se.aux))
}

## Non-Gaussian regression with fixed smoothing parameters
ngreg <- function(dc,family,sr,q,y,wt,offset,nu,alpha)
{
    nobs <- nrow(sr)
    nn <- ncol(sr)
    nxi <- nrow(q)
    nnull <- nn - nxi
    ## initialization
    cc <- dc[nnull+(1:nxi)]
    eta <- sr%*%dc
    if (!is.null(offset)) eta <- eta + offset
    if ((family=="nbinomial")&is.vector(y)) y <- cbind(y,nu[[1]])
    dev <- switch(family,
                  binomial=dev.resid.binomial(y,eta,wt),
                  nbinomial=dev.resid.nbinomial(y,eta,wt),
                  poisson=dev.resid.poisson(y,eta,wt),
                  Gamma=dev.resid.Gamma(y,eta,wt),
                  inverse.gaussian=dev.resid.inverse.gaussian(y,eta,wt),
                  weibull=dev.resid.weibull(y,eta,wt,nu[[1]]),
                  lognorm=dev0.resid.lognorm(y,eta,wt,nu[[1]]),
                  loglogis=dev0.resid.loglogis(y,eta,wt,nu[[1]]))
    dev <- sum(dev) + t(cc)%*%q%*%cc
    ## Newton iteration
    dc.new <- eta.new <- NULL
    dev.line <- function(x) {
        assign("dc.new",dc+x*dc.diff,inherits=TRUE)
        cc <- dc.new[nnull+(1:nxi)]
        eta.wk <- as.vector(sr%*%dc.new)
        if (!is.null(offset)) eta.wk <- eta.wk + offset
        assign("eta.new",eta.wk,inherits=TRUE)
        dev.wk <- switch(family,
                         binomial=dev.resid.binomial(y,eta.new,wt),
                         nbinomial=dev.resid.nbinomial(y,eta.new,wt),
                         poisson=dev.resid.poisson(y,eta.new,wt),
                         Gamma=dev.resid.Gamma(y,eta.new,wt),
                         inverse.gaussian=dev.resid.inverse.gaussian(y,eta.new,wt),
                         weibull=dev.resid.weibull(y,eta.new,wt,nu[[1]]),
                         lognorm=dev0.resid.lognorm(y,eta.new,wt,nu[[1]]),
                         loglogis=dev0.resid.loglogis(y,eta.new,wt,nu[[1]]))
        sum(dev.wk) + t(cc)%*%q%*%cc
    }
    iter <- 0
    flag <- 0
    flag2 <- 0
    repeat {
        iter <- iter+1
        dat <- switch(family,
                      binomial=mkdata.binomial(y,eta,wt,offset),
                      nbinomial=mkdata.nbinomial(y,eta,wt,offset,nu),
                      poisson=mkdata.poisson(y,eta,wt,offset),
                      Gamma=mkdata.Gamma(y,eta,wt,offset),
                      inverse.gaussian=mkdata.inverse.gaussian(y,eta,wt,offset),
                      weibull=mkdata.weibull(y,eta,wt,offset,nu),
                      lognorm=mkdata.lognorm(y,eta,wt,offset,nu),
                      loglogis=mkdata.loglogis(y,eta,wt,offset,nu))
        ## weighted least squares fit
        mumax <- 2*max(abs(t(sr)%*%dat$u+c(rep(0,nnull),q%*%dc[nnull+(1:nxi)])))
        w <- as.vector(sqrt(dat$wt))
        ywk <- w*dat$ywk
        srwk <- w*sr
        if (!is.finite(sum(w,ywk,srwk))) {
            if (flag) stop("gss error in gssanova: Newton iteration diverges")
            dc <- rep(0,nn)
            eta <- rep(0,nobs)
            if (!is.null(offset)) eta <- eta + offset
            dev <- switch(family,
                          binomial=dev.resid.binomial(y,eta,wt),
                          nbinomial=dev.resid.nbinomial(y,eta,wt),
                          poisson=dev.resid.poisson(y,eta,wt),
                          Gamma=dev.resid.Gamma(y,eta,wt),
                          inverse.gaussian=dev.resid.inverse.gaussian(y,eta,wt),
                          weibull=dev.resid.weibull(y,eta,wt,nu[[1]]),
                          lognorm=dev0.resid.lognorm(y,eta,wt,nu[[1]]),
                          loglogis=dev0.resid.loglogis(y,eta,wt,nu[[1]]))
            dev <- sum(dev)
            iter <- 0
            flag <- 1
            next
        }
        z <- .Fortran("reg",
                      as.double(srwk), as.integer(nobs), as.integer(nnull),
                      as.double(q), as.integer(nxi), as.double(ywk),
                      as.integer(4),
                      double(1), double(1), double(1), dc=double(nn),
                      as.double(.Machine$double.eps),
                      double(nn*nn), double(nn),
                      as.integer(c(rep(1,nnull),rep(0,nxi))),
                      double(max(nobs,nn)), integer(1), integer(1),
                      PACKAGE="gss")["dc"]
        dc.diff <- z$dc-dc
        repeat {
            dev.new <- dev.line(1)
            if (!is.finite(dev.new)) {
                dc.diff <- dc.diff/2
                next
            }
            if (!flag2) {
                if (dev.new-dev<1e-7*(1+abs(dev))) break
            }
            zz <- nlm0(dev.line,c(0,1),1e-3)
            dev.new <- dev.line(zz$est)
            break
        }
        disc0 <- max((mumax/(1+eta))^2,abs(eta.new-eta)/(1+eta))
        disc <- sum(dat$wt*((eta-eta.new)/(1+abs(eta)))^2)/sum(dat$wt)
        if (!is.finite(disc)) {
            if (flag) stop("gss error in gssanova: Newton iteration diverges")
            dc <- rep(0,nn)
            eta <- rep(0,nobs)
            if (!is.null(offset)) eta <- eta + offset
            dev <- switch(family,
                          binomial=dev.resid.binomial(y,eta,wt),
                          nbinomial=dev.resid.nbinomial(y,eta,wt),
                          poisson=dev.resid.poisson(y,eta,wt),
                          Gamma=dev.resid.Gamma(y,eta,wt),
                          inverse.gaussian=dev.resid.inverse.gaussian(y,eta,wt),
                          weibull=dev.resid.weibull(y,eta,wt,nu[[1]]),
                          lognorm=dev0.resid.lognorm(y,eta,wt,nu[[1]]),
                          loglogis=dev0.resid.loglogis(y,eta,wt,nu[[1]]))
            dev <- sum(dev)
            iter <- 0
            flag <- 1
            next
        }
        dc <- dc.new
        eta <- eta.new
        dev <- dev.new
        if (min(disc0,disc)<1e-7) break
        if (iter<=30) next
        if (!flag2) {
            flag2 <- 1
            iter <- 0
            next
        }
        warning("gss warning in gssanova: Newton iteration fails to converge")
        break
    }
    ## calculate cv
    dat <- switch(family,
                  binomial=mkdata.binomial(y,eta,wt,offset),
                  nbinomial=mkdata.nbinomial(y,eta,wt,offset,nu),
                  poisson=mkdata.poisson(y,eta,wt,offset),
                  Gamma=mkdata.Gamma(y,eta,wt,offset),
                  inverse.gaussian=mkdata.inverse.gaussian(y,eta,wt,offset),
                  weibull=mkdata.weibull(y,eta,wt,offset,nu),
                  lognorm=mkdata.lognorm(y,eta,wt,offset,nu),
                  loglogis=mkdata.loglogis(y,eta,wt,offset,nu))
    ## weighted least squares fit
    w <- as.vector(sqrt(dat$wt))
    ywk <- w*dat$ywk
    srwk <- w*sr
    z <- .Fortran("reg",
                  as.double(srwk), as.integer(nobs), as.integer(nnull),
                  as.double(q), as.integer(nxi), as.double(ywk),
                  as.integer(5),
                  double(1), double(1), double(1), dc=double(nn),
                  as.double(.Machine$double.eps),
                  chol=double(nn*nn), double(nn),
                  jpvt=as.integer(c(rep(1,nnull),rep(0,nxi))),
                  hat=double(max(nobs+1,nn)), rkv=integer(1), integer(1),
                  PACKAGE="gss")[c("dc","chol","jpvt","hat","rkv")]
    cv <- switch(family,
                 binomial=cv.binomial(y,eta,wt,z$hat[1:nobs],alpha),
                 poisson=cv.poisson(y,eta,wt,z$hat[1:nobs],alpha,sr,q),
                 Gamma=cv.Gamma(y,eta,wt,z$hat[1:nobs],z$hat[nobs+1],alpha),
                 inverse.gaussian=cv.inverse.gaussian(y,eta,wt,z$hat[1:nobs],z$hat[nobs+1],alpha),
                 nbinomial=cv.nbinomial(y,eta,wt,z$hat[1:nobs],alpha),
                 weibull=cv.weibull(y,eta,wt,z$hat[1:nobs],nu[[1]],alpha),
                 lognorm=cv.lognorm(y,eta,wt,z$hat[1:nobs],nu[[1]],alpha),
                 loglogis=cv.loglogis(y,eta,wt,z$hat[1:nobs],nu[[1]],alpha))
    c(z,cv,list(eta=eta))
}
