## Calculate Kullback-Leibler projection from gssanova objects
project.gssanova <- function(object,include,...)
{
    if (class(object)[1]=="gssanova0")
        stop("gss error: Kullback-Leibler projection is not implemented for gssanova0")
    nobs <- nrow(object$mf)
    nxi <- length(object$id.basis)
    labels.p <- object$lab.p
    ## evaluate full model
    family <- object$family
    eta <- object$eta
    y <- model.response(object$mf,"numeric")
    wt <- model.weights(object$mf)
    if(is.null(wt)) wt <- rep(1,nobs)
    offset <- model.offset(object$mf)
    if (!is.null(object$random)) {
        if (is.null(offset)) offset <- 0
        offset <- offset + object$random$z%*%object$b
    }
    nu <- object$nu
    y0 <- switch(family,
                 binomial=y0.binomial(y,eta,wt),
                 poisson=y0.poisson(eta),
                 Gamma=y0.Gamma(eta),
                 inverse.gaussian=y0.inverse.gaussian(eta),
                 nbinomial=y0.nbinomial(y,eta,nu),
                 weibull=y0.weibull(y,eta,nu),
                 lognorm=y0.lognorm(y,eta,nu),
                 loglogis=y0.loglogis(y,eta,nu))
    # calculate constant fit
    cfit <- switch(family,
                   binomial=cfit.binomial(y,wt,offset),
                   poisson=cfit.poisson(y,wt,offset),
                   Gamma=cfit.Gamma(y,wt,offset),
                   inverse.gaussian=cfit.inverse.gaussian(y,wt,offset),
                   nbinomial=cfit.nbinomial(y,wt,offset,nu),
                   weibull=cfit.weibull(y,wt,offset,nu),
                   lognorm=cfit.lognorm(y,wt,offset,nu),
                   loglogis=cfit.loglogis(y,wt,offset,nu))
    # calculate total entropy
    kl0 <- switch(family,
                  binomial=kl.binomial(eta,cfit,y0$wt),
                  poisson=kl.poisson(eta,cfit,wt),
                  Gamma=kl.Gamma(eta,cfit,wt),
                  inverse.gaussian=kl.inverse.gaussian(eta,cfit,wt),
                  nbinomial=kl.nbinomial(eta,cfit,wt,y0$nu),
                  weibull=kl.weibull(eta,cfit,wt,nu,y0$int),
                  lognorm=kl.lognorm(eta,cfit,wt,nu,y0),
                  loglogis=kl.loglogis(eta,cfit,wt,nu,y0))
    ## extract terms in subspace
    s <- matrix(1,nobs,1)
    philist <- object$term[["1"]]$iphi
    r <- NULL
    theta <- NULL
    nq.wk <- nq <- 0
    for (label in object$terms$labels) {
        if (label=="1") next
        if (label%in%labels.p) next
        x <- object$mf[,object$term[[label]]$vlist]
        x.basis <- object$mf[object$id.basis,object$term[[label]]$vlist]
        nphi <- object$term[[label]]$nphi
        nrk <- object$term[[label]]$nrk
        if (nphi) {
            phi <- object$term[[label]]$phi
            for (i in 1:nphi) {
                if (!any(label==include)) next
                philist <- c(philist,object$term[[label]]$iphi+(i-1))
                s <- cbind(s,phi$fun(x,nu=i,env=phi$env))
            }
        }
        if (nrk) {
            rk <- object$term[[label]]$rk
            for (i in 1:nrk) {
                nq.wk <- nq.wk + 1
                if (!any(label==include)) next
                nq <- nq + 1
                theta <- c(theta,object$theta[nq.wk])
                r <- array(c(r,rk$fun(x,x.basis,nu=i,env=rk$env,out=TRUE)),
                           c(nobs,nxi,nq))
            }
        }
    }
    if (!is.null(object$partial)) {
        nu <- length(object$d)-length(object$lab.p)
        matx.p <- model.matrix(object$partial$mt,object$mf)[,-1,drop=FALSE]
        matx.p <- scale(matx.p)
        for (label in labels.p) {
            nu <- nu+1
            if (!any(label==include)) next
            philist <- c(philist,nu)
            s <- cbind(s,matx.p[,label])
        }
    }
    ## calculate projection
    my.wls <- function(theta1=NULL) {
        if (!nq) {
            q <- matrix(0)
            sr <- cbind(s,0)
            z <- ngreg.proj(dc,family,sr,q,y0,wt,offset,nu)
        }
        else {
            theta.wk <- 1:nq
            theta.wk[fix] <- theta[fix]
            if (nq-1) theta.wk[-fix] <- theta1
            sr <- 0
            for (i in 1:nq) sr <- sr + 10^theta.wk[i]*r[,,i]
            q <- sr[object$id.basis,]
            sr <- cbind(s,sr)
            z <- ngreg.proj(dc,family,sr,q,y0,wt,offset,nu)
        }
        assign("dc",z$dc,inherits=TRUE)
        assign("eta1",z$eta,inherits=TRUE)
        z$kl
    }
    cv.wk <- function(theta) cv.scale*my.wls(theta)+cv.shift
    ## initialization
    if (nq) {
        r.wk <- 0
        for (i in 1:nq) r.wk <- r.wk + 10^theta[i]*r[,,i]
        if (is.null(s)) theta.wk <- 0
        else theta.wk <- log10(sum(s^2)/ncol(s)/sum(r.wk^2)*nxi) / 2
        theta <- theta + theta.wk
        tmp <- NULL
        for (i in 1:nq) tmp <- c(tmp,10^theta[i]*sum(r[cbind(object$id.basis,1:nxi,i)]))
        fix <- rev(order(tmp))[1]
    }
    ## projection
    if (nq) dc <- c(object$d[philist],10^(-theta.wk)*object$c)
    else dc <- c(object$d[philist],0)
    eta1 <- NULL
    if (nq>1) {
        if (object$skip.iter) kl <- my.wls(theta[-fix])
        else {
            ## scale and shift cv
            tmp <- abs(my.wls(theta[-fix]))
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
            zz <- nlm(cv.wk,theta[-fix],stepmax=1,ndigit=7)
            kl <- my.wls(zz$est)
        }
    }
    else kl <- my.wls()
    ## check
    kl1 <- switch(family,
                  binomial=kl.binomial(eta1,cfit,y0$wt),
                  poisson=kl.poisson(eta1,cfit,wt),
                  Gamma=kl.Gamma(eta1,cfit,wt),
                  inverse.gaussian=kl.inverse.gaussian(eta1,cfit,wt),
                  nbinomial=kl.nbinomial(eta1,cfit,wt,y0$nu),
                  weibull=kl.weibull(eta1,cfit,wt,nu,y0$int),
                  lognorm=kl.lognorm(eta1,cfit,wt,nu,y0),
                  loglogis=kl.loglogis(eta1,cfit,wt,nu,y0))
    list(ratio=kl/kl0,kl=kl,check=(kl+kl1)/kl0)
}

## KL projection with Non-Gaussian regression
ngreg.proj <- function(dc,family,sr,q,y0,wt,offset,nu)
{
    ## initialization
    q <- 10^(-5)*q
    eta <- sr%*%dc
    nobs <- length(eta)
    nn <- ncol(as.matrix(sr))
    nxi <- ncol(q)
    nnull <- nn-nxi
    if (!is.null(offset)) eta <- eta + offset
    fit1 <- switch(family,
                   binomial=proj0.binomial(y0,eta,offset),
                   poisson=proj0.poisson(y0,eta,wt,offset),
                   Gamma=proj0.Gamma(y0,eta,wt,offset),
                   inverse.gaussian=proj0.inverse.gaussian(y0,eta,wt,offset),
                   nbinomial=proj0.nbinomial(y0,eta,wt,offset),
                   weibull=proj0.weibull(y0,eta,wt,offset,nu),
                   lognorm=proj0.lognorm(y0,eta,wt,offset,nu),
                   loglogis=proj0.loglogis(y0,eta,wt,offset,nu))
    kl <- fit1$kl+t(dc[-(1:nnull)])%*%q%*%dc[-(1:nnull)]/2
    ## Newton iteration
    dc.new <- eta.new <- NULL
    kl.line <- function(x) {
        assign("dc.new",dc+x*dc.diff,inherits=TRUE)
        eta.wk <- sr%*%dc.new
        if (!is.null(offset)) eta.wk <- eta.wk + offset
        assign("eta.new",eta.wk,inherits=TRUE)
        fit.wk <- switch(family,
                         binomial=proj0.binomial(y0,eta.new,offset),
                         poisson=proj0.poisson(y0,eta.new,wt,offset),
                         Gamma=proj0.Gamma(y0,eta.new,wt,offset),
                         inverse.gaussian=proj0.inverse.gaussian(y0,eta.new,wt,offset),
                         nbinomial=proj0.nbinomial(y0,eta.new,wt,offset),
                         weibull=proj0.weibull(y0,eta.new,wt,offset,nu),
                         lognorm=proj0.lognorm(y0,eta.new,wt,offset,nu),
                         loglogis=proj0.loglogis(y0,eta.new,wt,offset,nu))
        assign("fit1",fit.wk,inherits=TRUE)
        fit1$kl+t(dc.new[-(1:nnull)])%*%q%*%dc.new[-(1:nnull)]/2
    }
    iter <- 0
    flag <- 0
    flag2 <- 0
    repeat {
        iter <- iter+1
        ## weighted least squares fit
        if (!is.finite(sum(fit1$wt,fit1$ywk))) {
            if (flag) stop("gss error in project.gssanova: Newton iteration diverges")
            dc <- rep(0,nn)
            eta <- rep(0,nobs)
            if (!is.null(offset)) eta <- eta + offset
            fit1 <- switch(family,
                           binomial=proj0.binomial(y0,eta,offset),
                           poisson=proj0.poisson(y0,eta,wt,offset),
                           Gamma=proj0.Gamma(y0,eta,wt,offset),
                           inverse.gaussian=proj0.inverse.gaussian(y0,eta,wt,offset),
                           nbinomial=proj0.nbinomial(y0,eta,wt,offset),
                           weibull=proj0.weibull(y0,eta,wt,offset,nu),
                           lognorm=proj0.lognorm(y0,eta,wt,offset,nu),
                           loglogis=proj0.loglogis(y0,eta,wt,offset,nu))
            kl <- fit1$kl
            iter <- 0
            flag <- 1
            next
        }
        mumax <- max(abs(t(sr)%*%fit1$u+c(rep(0,nnull),q%*%dc[-(1:nnull)])))
        w <- sqrt(as.vector(fit1$wt))
        z <- .Fortran("reg",
                      as.double(w*sr), as.integer(nobs), as.integer(nnull),
                      as.double(q), as.integer(nxi), as.double(w*fit1$ywk),
                      as.integer(4),
                      double(1), double(1), double(1), dc=double(nn),
                      as.double(.Machine$double.eps),
                      double(nn*nn), double(nn), as.integer(rep(0,nn)),
                      double(max(nobs,nn)), integer(1), integer(1),
                      PACKAGE="gss")["dc"]
        dc.diff <- z$dc-dc
        repeat {
            kl.new <- kl.line(1)
            if (!is.finite(kl.new)) {
                dc.diff <- dc.diff/2
                next
            }
            if (!flag2) {
                if (kl.new-kl<1e-7*(1+abs(kl))) break
            }
            zz <- nlm0(kl.line,c(0,1),1e-3)
            kl.new <- kl.line(zz$est)
            break
        }
        disc0 <- max((mumax/(1+kl))^2,abs(kl.new-kl)/(1+kl))
        disc <- sum(fit1$wt*((eta-eta.new)/(1+abs(eta)))^2)/sum(fit1$wt)
        if (is.nan(disc)) {
            if (flag) stop("gss error in project.gssanova: Newton iteration diverges")
            dc <- rep(0,nn)
            eta <- rep(0,nobs)
            if (!is.null(offset)) eta <- eta + offset
            fit1 <- switch(family,
                           binomial=proj0.binomial(y0,eta,offset),
                           poisson=proj0.poisson(y0,eta,wt,offset),
                           Gamma=proj0.Gamma(y0,eta,wt,offset),
                           inverse.gaussian=proj0.inverse.gaussian(y0,eta,wt,offset),
                           nbinomial=proj0.nbinomial(y0,eta,wt,offset),
                           weibull=proj0.weibull(y0,eta,wt,offset,nu),
                           lognorm=proj0.lognorm(y0,eta,wt,offset,nu),
                           loglogis=proj0.loglogis(y0,eta,wt,offset,nu))
            kl <- fit1$kl
            iter <- 0
            flag <- 1
            next
        }
        dc <- dc.new
        eta <- eta.new
        kl <- kl.new
        if (min(disc0,disc)<1e-5) break
        if (iter<=30) next
        if (!flag2) {
            flag2 <- 1
            iter <- 0
            next
        }
        warning("gss warning in gssanova: Newton iteration fails to converge")
        break
    }
    fit1 <- switch(family,
                   binomial=proj0.binomial(y0,eta,offset),
                   poisson=proj0.poisson(y0,eta,wt,offset),
                   Gamma=proj0.Gamma(y0,eta,wt,offset),
                   inverse.gaussian=proj0.inverse.gaussian(y0,eta,wt,offset),
                   nbinomial=proj0.nbinomial(y0,eta,wt,offset),
                   weibull=proj0.weibull(y0,eta,wt,offset,nu),
                   lognorm=proj0.lognorm(y0,eta,wt,offset,nu),
                   loglogis=proj0.loglogis(y0,eta,wt,offset,nu))
    kl <- fit1$kl
    list(dc=dc,eta=eta,kl=kl)
}
