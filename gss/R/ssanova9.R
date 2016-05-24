## Fit ssanova model with correlated data
ssanova9 <- function(formula,type=NULL,data=list(),subset,
                     offset,na.action=na.omit,partial=NULL,
                     method="v",alpha=1.4,varht=1,
                     id.basis=NULL,nbasis=NULL,seed=NULL,cov,
                     skip.iter=FALSE)
{
    ## Obtain model frame and model terms
    mf <- match.call()
    mf$type <- mf$method <- mf$varht <- mf$partial <- NULL
    mf$alpha <- mf$id.basis <- mf$nbasis <- mf$seed <- NULL
    mf$cov <- mf$skip.iter <- NULL
    mf[[1]] <- as.name("model.frame")
    mf <- eval(mf,parent.frame())
    ## Generate sub-basis
    nobs <- dim(mf)[1]
    if (is.null(id.basis)) {
        if (is.null(nbasis))  nbasis <- max(30,ceiling(10*nobs^(2/9)))
        if (nbasis>=nobs)  nbasis <- nobs
        if (!is.null(seed))  set.seed(seed)
        id.basis <- sample(nobs,nbasis)
    }
    else {
        if (max(id.basis)>nobs|min(id.basis)<1)
            stop("gss error in ssanova9: id.basis out of range")
        nbasis <- length(id.basis)
    }
    ## Generate terms
    term <- mkterm(mf,type)
    ## Generate cov
    if (is.null(cov$fun)) {
        type <- cov[[1]]
        if (type=="arma") {
            pq <- cov[[2]]
            cov <- mkcov.arma(pq[1],pq[2],nobs)
        }
        if (type=="long") {
            if (nobs<length(cov[[2]])) id <- cov[[2]][subset]
            else id <- cov[[2]]
            cov <- mkcov.long(id)
        }
        if (type=="known") {
            cov <- mkcov.known(as.matrix(cov[[2]]))
        }
    }        
    ## Generate s and r
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
        stop("gss error in ssanova9: use lm for models with only unpenalized terms")
    ## Add the partial term
    if (!is.null(partial)) {
        mf.p <- model.frame(partial,data)
        for (lab in colnames(mf.p)) mf[,lab] <- mf.p[,lab]
        mt.p <- attr(mf.p,"terms")
        lab.p <- labels(mt.p)
        matx.p <- model.matrix(mt.p,data)[,-1,drop=FALSE]
        if (dim(matx.p)[1]!=dim(mf)[1])
            stop("gss error in ssanova9: partial data are of wrong size")
        matx.p <- scale(matx.p)
        center.p <- attr(matx.p,"scaled:center")
        scale.p <- attr(matx.p,"scaled:scale")
        s <- cbind(s,matx.p)
        part <- list(mt=mt.p,center=center.p,scale=scale.p)
    }
    else part <- lab.p <- NULL
    if (qr(s)$rank<dim(s)[2])
        stop("gss error in ssanova9: unpenalized terms are linearly dependent")
    ## Prepare the data
    y <- model.response(mf,"numeric")
    offset <- model.offset(mf)
    if (!is.null(offset)) {
        term$labels <- c(term$labels,"offset")
        term$offset <- list(nphi=0,nrk=0)
        y <- y - offset
    }
    ## Fit the model
    if (nq==1) {
        r <- r[,,1]
        z <- sspreg91(s,r,r[id.basis,],y,cov,method,alpha,varht)
    }
    else z <- mspreg91(s,r,id.basis,y,cov,method,alpha,varht,skip.iter)
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
    obj <- c(list(call=match.call(),mf=mf,terms=term,desc=desc,alpha=alpha,
                  id.basis=id.basis,partial=part,lab.p=lab.p,cov=cov,
                  skip.iter=skip.iter),z)
    class(obj) <- c("ssanova9","ssanova")
    obj
}

## Fit Single Smoothing Parameter (Gaussian) REGression
sspreg91 <- function(s,r,q,y,cov,method,alpha,varht)
{
    ## get dimensions
    nobs <- nrow(r)
    nxi <- ncol(r)
    if (!is.null(s)) {
        if (is.vector(s)) nnull <- 1
        else nnull <- ncol(s)
    }
    else nnull <- 0
    nn <- nxi + nnull
    ## cv function
    cv <- function(lambda) {
        q.wk <- 10^(lambda[1]+theta)*q
        if (length(lambda)-1) {
            ww <- cov$fun(lambda[-1],cov$env)
            ww <- chol(ww)
            y.wk <- forwardsolve(t(ww),y)
            s.wk <- forwardsolve(t(ww),s)
            r.wk <- forwardsolve(t(ww),r)
        }
        z <- .Fortran("reg",
                      as.double(cbind(s.wk,10^theta*r.wk)),
                      as.integer(nobs), as.integer(nnull),
                      as.double(q.wk), as.integer(nxi), as.double(y.wk),
                      as.integer(switch(method,"u"=1,"v"=2,"m"=3)),
                      as.double(alpha), varht=as.double(varht),
                      score=double(1), dc=double(nn),
                      as.double(.Machine$double.eps),
                      chol=double(nn*nn), double(nn),
                      jpvt=as.integer(c(rep(1,nnull),rep(0,nxi))),
                      wk=double(3*nobs+nnull), rkv=integer(1), info=integer(1),
                      PACKAGE="gss")[c("score","varht","dc","chol","jpvt","wk","rkv","info")]
        if (z$info) stop("gss error in ssanova9: evaluation of GML score fails")
        if (!nnull|method%in%c("u","v")) detw <- 2*sum(log(diag(ww)))
        else {
            wk <- qr.qty(qr(s),t(ww))[-(1:nnull),]
            detw <- sum(log(eigen(wk%*%t(wk))$value))
        }
        if (method=="m") score <- z$score*exp(detw/(nobs-nnull))
        if (method=="u") score <- z$wk[1]/varht+detw/nobs+2*alpha*z$wk[2]
        if (method=="v") score <- log(z$wk[1])+detw/nobs+2*alpha*z$wk[2]/(1-z$wk[2])
        alpha.wk <- max(0,log.la0-lambda[1]-5)*(3-alpha) + alpha
        alpha.wk <- min(alpha.wk,3)
        if (alpha.wk>alpha) {
            if (method=="u") score <- score + (alpha.wk-alpha)*2*z$wk[2]
            if (method=="v") score <- score + (alpha.wk-alpha)*2*z$wk[2]/(1-z$wk[2])
        }
        z$score <- score
        assign("fit",z[c(1:5,7)],inherits=TRUE)
        score
    }
    cv.wk <- function(lambda) cv.scale*cv(lambda)+cv.shift
    ## initialization
    tmp <- sum(r^2)
    if (is.null(s)) theta <- 0
    else theta <- log10(sum(s^2)/nnull/tmp*nxi) / 2
    log.la0 <- log10(tmp/sum(diag(q))) + theta
    ## lambda search
    fit <- NULL
    la <- c(log.la0,cov$init)
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
                warning("gss warning in ssanova9: iteration for model selection fails to converge")
                break
            }
        }
    }
    else {
        ww <- cov$fun(cov$env)
        ww <- chol(ww)
        y.wk <- forwardsolve(t(ww),y)
        s.wk <- forwardsolve(t(ww),s)
        r.wk <- forwardsolve(t(ww),r)
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
    lambda <- zz$est
    jk1 <- cv(lambda)
    q.wk <- 10^(theta)*q
    if (length(lambda)-1) ww <- cov$fun(lambda[-1],cov$env)
    else ww <- cov$fun(cov$env)
    ww <- chol(ww)
    s.wk <- forwardsolve(t(ww),s)
    r.wk <- forwardsolve(t(ww),r)
    se.aux <- regaux(s.wk,10^theta*r.wk,q.wk,lambda[1],fit)
    c <- fit$dc[nnull+(1:nxi)]
    if (nnull) d <- fit$dc[1:nnull]
    else d <- NULL
    c(list(method=method,theta=theta,c=c,d=d,nlambda=lambda[1],zeta=lambda[-1]),
      fit[-3],list(se.aux=se.aux))
}

## Fit Multiple Smoothing Parameter (Gaussian) REGression
mspreg91 <- function(s,r,id.basis,y,cov,method,alpha,varht,skip.iter)
{
    ## get dimensions
    nobs <- nrow(r)
    nxi <- ncol(r)
    if (!is.null(s)) {
        if (is.vector(s)) nnull <- 1
        else nnull <- ncol(s)
    }
    else nnull <- 0
    nn <- nxi + nnull
    nq <- dim(r)[3]
    ## cv function
    cv <- function(theta) {
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
        q.wk <- 10^nlambda*r.wk0[id.basis,]
        if (length(theta)-nq) {
            ww <- cov$fun(theta[-(1:nq)],cov$env)
            ww <- chol(ww)
            y.wk <- forwardsolve(t(ww),y)
            s.wk <- forwardsolve(t(ww),s)
        }
        r.wk0 <- forwardsolve(t(ww),r.wk0)
        z <- .Fortran("reg",
                      as.double(cbind(s.wk,r.wk0)), as.integer(nobs), as.integer(nnull),
                      as.double(q.wk), as.integer(nxi), as.double(y.wk),
                      as.integer(switch(method,"u"=1,"v"=2,"m"=3)),
                      as.double(alpha), varht=as.double(varht),
                      score=double(1), dc=double(nn),
                      as.double(.Machine$double.eps),
                      chol=double(nn*nn), double(nn),
                      jpvt=as.integer(c(rep(1,nnull),rep(0,nxi))),
                      wk=double(3*nobs+nnull), rkv=integer(1), info=integer(1),
                      PACKAGE="gss")[c("score","varht","dc","chol","jpvt","wk","rkv","info")]
        if (z$info) stop("gss error in ssanova9: evaluation of GML score fails")
        if (!nnull|method%in%c("u","v")) detw <- 2*sum(log(diag(ww)))
        else {
            wk <- qr.qty(qr(s),t(ww))[-(1:nnull),]
            detw <- sum(log(eigen(wk%*%t(wk))$value))
        }
        if (method=="m") score <- z$score*exp(detw/(nobs-nnull))
        if (method=="u") score <- z$wk[1]/varht+detw/nobs+2*alpha*z$wk[2]
        if (method=="v") score <- log(z$wk[1])+detw/nobs+2*alpha*z$wk[2]/(1-z$wk[2])
        alpha.wk <- max(0,theta[1:nq]-log.th0-5)*(3-alpha) + alpha
        alpha.wk <- min(alpha.wk,3)
        if (alpha.wk>alpha) {
            if (method=="u") score <- score + (alpha.wk-alpha)*2*z$wk[2]
            if (method=="v") score <- score + (alpha.wk-alpha)*2*z$wk[2]/(1-z$wk[2])
        }
        z$score <- score
        assign("fit",z[c(1:5,7)],inherits=TRUE)
        score
    }
    cv.wk <- function(theta) cv.scale*cv(theta)+cv.shift
    ## initialization
    theta <- -log10(apply(r[id.basis,,],3,function(x)sum(diag(x))))
    r.wk <- 0
    for (i in 1:nq) {
        r.wk <- r.wk + 10^theta[i]*r[,,i]
    }
    ## theta adjustment
    z <- sspreg91(s,r.wk,r.wk[id.basis,],y,cov,method,alpha,varht)
    theta <- theta + z$theta
    r.wk <- 0
    for (i in 1:nq) {
        theta[i] <- 2*theta[i] + log10(t(z$c)%*%r[id.basis,,i]%*%z$c)
        r.wk <- r.wk + 10^theta[i]*r[,,i]
    }
    log.la0 <- log10(sum(r.wk^2)/sum(diag(r.wk[id.basis,])))
    log.th0 <- theta-log.la0
    ## lambda search
    z <- sspreg91(s,r.wk,r.wk[id.basis,],y,cov,method,alpha,varht)
    nlambda <- z$nlambda
    log.th0 <- log.th0 + z$nlambda
    theta <- theta + z$theta
    ## early return
    if (skip.iter) {
        z$theta <- theta
        return(z)
    }
    ## theta search
    fit <- NULL
    counter <- 0
    r.wk <- 0
    for (i in 1:nq) {
        r.wk <- r.wk + 10^theta[i]*r[,,i]
    }
    theta.old <- theta
    theta <- c(theta,z$zeta)
    ## scale and shift cv
    if (length(theta)==nq) {
        ww <- cov$fun(cov$env)
        ww <- chol(ww)
        y.wk <- forwardsolve(t(ww),y)
        s.wk <- forwardsolve(t(ww),s)
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
            warning("gss warning in ssanova9: iteration for model selection fails to converge")
            break
        }
    }
    ## return
    theta <- zz$est
    jk1 <- cv(theta)
    r.wk <- 0
    for (i in 1:nq) {
        r.wk <- r.wk + 10^theta[i]*r[,,i]
    }
    q.wk <- r.wk[id.basis,]
    if (length(theta)-nq) ww <- cov$fun(theta[-(1:nq)],cov$env)
    else ww <- cov$fun(cov$env)
    ww <- chol(ww)
    s.wk <- forwardsolve(t(ww),s)
    r.wk <- forwardsolve(t(ww),r.wk)
    se.aux <- regaux(s.wk,r.wk,q.wk,nlambda,fit)
    c <- fit$dc[nnull+(1:nxi)]
    if (nnull) d <- fit$dc[1:nnull]
    else d <- NULL
    c(list(method=method,theta=theta[1:nq],c=c,d=d,nlambda=nlambda,
           zeta=theta[-(1:nq)]),fit[-3],list(se.aux=se.aux))
}

mkcov.arma <- function(p,q,n) {
    ## check inputs
    if ((p<0)|(q<0))
      stop("gss error in mkcov.arma: ARMA orders must be non-negative")
    if (p+q==0) stop("gss error in mkcov.arma: use ssanova for independent data")
    if (n<(p+1)) stop("gss error in mkcov.arma: AR order too high")
    env <- list(p=p,q=q,n=n)
    init <- rep(0,p+q)
    fun <- function(x,env) {
        p <- env$p
        q <- env$q
        n <- env$n
        ## arma coefficients
        if (p) {
            a <- (1-exp(-x[1:p]))/(1+exp(-x[1:p]))
            if (p-1) {
                for (j in 2:p) a[1:(j-1)] <- a[1:(j-1)]-a[j]*a[(j-1):1]
            }
        }
        else a <- NULL
        if (q) {
            b <- (1-exp(-x[p+(1:q)]))/(1+exp(-x[p+(1:q)]))
            if (q-1) {
                for (j in 2:q) b[1:(j-1)] <- b[1:(j-1)]-b[j]*b[(j-1):1]
            }
        }
        else b <- NULL
        ## psi
        psi <- 1
        if (qq <- max(p-1,q)) {
            for(i in 1:qq) {
                wk <- ifelse(i<=q,-b[i],0)
                if(p) {
                    for(j in 1:min(i,p)) wk <- wk + a[j]*psi[i-j+1]
                }
                psi<-c(psi,wk)
            }
        }
        ## autocovariance
        aa <- bb <- 1
        if (p) aa <- c(aa,-a)
        if (q) bb <- c(bb,-b)
        if (length(bb)<p) bb <- c(bb,rep(0,p-length(bb)))
        wk <- rep(0,n)
        for (i in 0:qq) {
            wk[i+1] <- sum(bb[(i:qq)+1]*psi[1:(qq-i+1)])
        }
        h <- matrix(0,p+1,p+1)
        for (i in 1:(p+1)) {
            for (j in 1:(p+1)) {
                jj <- abs(i-j)+1
                h[i,jj] <- h[i,jj]+aa[j]
            }
        }
        gg <- solve(h,wk[1:(p+1)])
        for (i in (p+2):n) {
            wk1 <- wk[i]
            if (p) wk1 <- wk1 + sum(a*rev(gg)[1:p])
            gg <- c(gg,wk1)
        }
        ## return
        cov <- matrix(0,n,n)
        for (i in 1:n) {
            wk1 <- diag(gg[i],n-i+1)
            cov[1:(n-i+1),i:n] <- cov[1:(n-i+1),i:n]+wk1
        }
        cov
    }
    list(fun=fun,env=env,init=init)
}
para.arma <- function(fit) {
    x <- fit$zeta
    p <- fit$cov$env$p
    q <- fit$cov$env$q
    n <- fit$cov$env$n
    ## arma coefficients
    if (p) {
        a <- (1-exp(-x[1:p]))/(1+exp(-x[1:p]))
        if (p-1) {
            for (j in 2:p) a[1:(j-1)] <- a[1:(j-1)]-a[j]*a[(j-1):1]
        }
    }
    else a <- NULL
    if (q) {
        b <- (1-exp(-x[p+(1:q)]))/(1+exp(-x[p+(1:q)]))
        if (q-1) {
            for (j in 2:q) b[1:(j-1)] <- b[1:(j-1)]-b[j]*b[(j-1):1]
        }
    }
    else b <- NULL
    list(a=a,b=b)
}

mkcov.long <- function(id) {
    ## check inputs
    if (!is.factor(id)) stop("gss error in mkcov.long: subject id must be a factor")
    n <- length(id)
    lvl <- levels(id)
    z <- NULL
    for (i in lvl) z <- cbind(z,as.numeric(id==i))
    env <- list(zz=z%*%t(z),n=n)
    init <- -log(n)
    fun <- function(x,env) {
        zz <- env$zz
        n <- env$n
        diag(1,n)+exp(x)*zz
    }
    list(fun=fun,env=env,init=init)
}

mkcov.known <- function(w) {
    ## check inputs
    if (nrow(w)-ncol(w)) stop("gss error in mkcov.known: matrix is not square")
    env <- w
    init <- NULL
    fun <- function(env) env
    list(fun=fun,env=env,init=init)
}
