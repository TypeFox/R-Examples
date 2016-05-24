## Calculate Kullback-Leibler projection from sscden objects
project.sscden <- function(object,include,...)
{
    mf <- object$mf
    term <- object$term
    id.basis <- object$id.basis
    qd.pt <- object$yquad$pt
    qd.wt <- object$yquad$wt
    xx.wt <- object$xx.wt
    ## evaluate full model
    x <- object$mf[!object$x.dup.ind,object$xnames,drop=FALSE]
    fit0 <- object$fit
    ## extract terms in subspace
    include <- union(object$ynames,include)
    nmesh <- length(qd.wt)
    nbasis <- length(id.basis)
    nx <- length(xx.wt)
    qd.s <- NULL
    qd.r <- as.list(NULL)
    theta <- d <- q <- NULL
    nu.wk <- nu <- nq.wk <- nq <- 0
    for (label in term$labels) {
        vlist <- term[[label]]$vlist
        x.list <- object$xnames[object$xnames%in%vlist]
        y.list <- object$ynames[object$ynames%in%vlist]
        xy.basis <- mf[id.basis,vlist]
        qd.xy <- data.frame(matrix(0,nmesh,length(vlist)))
        names(qd.xy) <- vlist
        qd.xy[,y.list] <- qd.pt[,y.list]
        if (length(x.list)) xx <- x[,x.list,drop=FALSE]
        else xx <- NULL
        nphi <- term[[label]]$nphi
        nrk <- term[[label]]$nrk
        if (nphi) {
            phi <- term[[label]]$phi
            for (i in 1:nphi) {
                nu.wk <- nu.wk+1
                if (is.null(xx)) {
                    if (!any(label==include)) next
                    nu <- nu+1
                    d <- c(d,object$d[nu.wk])
                    s.wk <- phi$fun(qd.xy[,,drop=TRUE],nu=i,env=phi$env)
                    wk <- matrix(s.wk,nmesh,nx)
                    qd.s <- array(c(qd.s,wk),c(nmesh,nx,nu))
                }
                else {
                    if (!any(label==include)) next
                    nu <- nu+1
                    d <- c(d,object$d[nu.wk])
                    wk <- NULL
                    for (j in 1:nx) {
                        qd.xy[,x.list] <- xx[rep(j,nmesh),]
                        wk <- cbind(wk,phi$fun(qd.xy,i,phi$env))
                    }
                    qd.s <- array(c(qd.s,wk),c(nmesh,nx,nu))
                }
            }
        }
        if (nrk) {
            rk <- term[[label]]$rk
            for (i in 1:nrk) {
                nq.wk <- nq.wk+1
                if (is.null(xx)) {
                    if (!any(label==include)) next
                    nq <- nq+1
                    theta <- c(theta,object$theta[nq.wk])
                    qd.r.wk <- rk$fun(qd.xy[,,drop=TRUE],xy.basis,nu=i,env=rk$env,out=TRUE)
                    qd.r[[nq]] <- qd.r.wk
                    q <- cbind(q,rk$fun(xy.basis,xy.basis,i,rk$env,out=FALSE))
                }
                else {
                    if (!any(label==include)) next
                    nq <- nq+1
                    theta <- c(theta,object$theta[nq.wk])
                    qd.wk <- NULL
                    for (j in 1:nx) {
                        qd.xy[,x.list] <- xx[rep(j,nmesh),]
                        qd.wk <- array(c(qd.wk,rk$fun(qd.xy,xy.basis,i,rk$env,TRUE)),
                                       c(nmesh,nbasis,j))
                    }
                    qd.r[[nq]] <- qd.wk
                    q <- cbind(q,rk$fun(xy.basis,xy.basis,i,rk$env,out=FALSE))
                }
            }
        }
    }
    if (is.null(qd.s)&is.null(qd.r))
        stop("gss error in project.sscden: include some terms")
    nnull <- length(d)
    nxis <- nbasis+nnull
    ## calculate projection
    rkl <- function(theta1=NULL) {
        theta.wk <- 1:nq
        theta.wk[fix] <- theta[fix]
        if (nq-1) theta.wk[-fix] <- theta1
        qd.rs <- array(0,c(nmesh,nbasis,nx))
        for (i in 1:nq) {
            if (length(dim(qd.r[[i]]))==3) qd.rs <- qd.rs + 10^theta[i]*qd.r[[i]]
            else qd.rs <- qd.rs + as.vector(10^theta[i]*qd.r[[i]])
        }
        qd.rs <- aperm(qd.rs,c(1,3,2))
        qd.rs <- array(c(qd.rs,qd.s),c(nmesh,nx,nxis))
        qd.rs <- aperm(qd.rs,c(1,3,2))
        z <- .Fortran("cdenrkl",
                      cd=as.double(cd), as.integer(nxis),
                      as.double(qd.rs), as.integer(nmesh), as.integer(nx),
                      as.double(xx.wt), as.double(qd.wt), as.double(t(fit0)),
                      as.double(.Machine$double.eps),
                      wt=double(nmesh*nx), double(nmesh*nx), double(nxis),
                      double(nxis), double(nxis*nxis), double(nxis*nxis),
                      integer(nxis), double(nxis), as.double(1e-6), as.integer(30),
                      info=integer(1), PACKAGE="gss")
        if (z$info==1)
            stop("gss error in project.sscden: Newton iteration diverges")
        if (z$info==2)
            warning("gss warning in project.sscden: Newton iteration fails to converge")
        assign("cd",z$cd,inherits=TRUE)
        z$wt[1]
    }
    cv.wk <- function(theta) cv.scale*rkl(theta)+cv.shift
    if (nq) {
        ## initialization
        if (!nnull) theta.wk <- 0
        else {
            qd.r.wk <- array(0,c(nmesh,nbasis,nx))
            for (i in 1:nq) {
                if (length(dim(qd.r[[i]]))==3) qd.r.wk <- qd.r.wk + 10^theta[i]*qd.r[[i]]
                else qd.r.wk <- qd.r.wk + as.vector(10^theta[i]*qd.r[[i]])
            }
            v.s <- v.r <- 0
            for (i in 1:nx) {
                mu.s <- apply(fit0[i,]*qd.s[,i,,drop=FALSE],2,sum)
                v.s.wk <- apply(fit0[i,]*qd.s[,i,,drop=FALSE]^2,2,sum)-mu.s^2
                mu.r <- apply(fit0[i,]*qd.r.wk[,,i,drop=FALSE],2,sum)
                v.r.wk <- apply(fit0[i,]*qd.r.wk[,,i,drop=FALSE]^2,2,sum)-mu.r^2
                v.s <- v.s + xx.wt[i]*v.s.wk
                v.r <- v.r + xx.wt[i]*v.r.wk
            }
            theta.wk <- log10(sum(v.s)/nnull/sum(v.r)*nbasis) / 2
        }
        theta <- theta + theta.wk
        tmp <- NULL
        for (i in 1:nq) tmp <- c(tmp,10^theta[i]*sum(q[,i]))
        fix <- rev(order(tmp))[1]
        ## projection
        cd <- c(10^(-theta.wk)*object$c,d)
        mesh1 <- NULL
        if (nq-1) {
            if (object$skip.iter) kl <- rkl(theta[-fix])
            else {
                if (nq-2) {
                    ## scale and shift cv
                    tmp <- abs(rkl(theta[-fix]))
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
                    zz <- nlm(cv.wk,theta[-fix],stepmax=.5,ndigit=7)
                }
                else {
                    the.wk <- theta[-fix]
                    repeat {
                        mn <- the.wk-1
                        mx <- the.wk+1
                        zz <- nlm0(rkl,c(mn,mx))
                        if (min(zz$est-mn,mx-zz$est)>=1e-3) break
                        else the.wk <- zz$est
                    }
                }
                kl <- rkl(zz$est)
            }
        }
        else kl <- rkl()
    }
    else {
        z <- .Fortran("cdenrkl",
                      cd=as.double(d), as.integer(nnull),
                      as.double(aperm(qd.s,c(1,3,2))), as.integer(nmesh), as.integer(nx),
                      as.double(xx.wt), as.double(qd.wt), as.double(t(fit0)),
                      as.double(.Machine$double.eps),
                      wt=double(nmesh*nx), double(nmesh*nx), double(nnull),
                      double(nnull), double(nnull*nnull), double(nnull*nnull),
                      integer(nnull), double(nnull), as.double(1e-6), as.integer(30),
                      info=integer(1), PACKAGE="gss")
        if (z$info==1)
            stop("gss error in project.sscden: Newton iteration diverges")
        if (z$info==2)
            warning("gss warning in project.sscden: Newton iteration fails to converge")
        kl <- z$wt[1]
    }
    ## cfit
    cfit <- matrix(1,nmesh,nx)
    for (ylab in object$ynames) {
        y <- object$mf[[ylab]]
        if (is.factor(y)) {
            lvl <- levels(y)
            if (is.null(object$cnt)) wk <- table(y)
            else wk <- table(rep(y,object$cnt))
            wk <- wk/sum(wk)
            nlvl <- length(wk)
            for (j in 1:nlvl) {
                id <- (1:nmesh)[qd.pt[,ylab]==lvl[j]]
                cfit[id,] <- cfit[id,]*wk[j]
            }
        }
        else {
            if (!is.vector(y)) qd.wk <- object$yquad
            else qd.wk <- NULL
            qd.wk <- object$yquad
            form <- as.formula(paste("~",ylab))
            wk <- ssden(form,data=object$mf,quad=qd.wk,
                        domain=object$ydomain,alpha=object$alpha,
                        id.basis=object$id.basis)
            cfit <- cfit*dssden(wk,qd.pt[ylab])
        }
    }
    cfit <- t(cfit*qd.wt)
    ## return
    kl0 <- 0
    for (i in 1:nx) {
        wk <- sum(log(fit0[i,]/cfit[i,])*fit0[i,])
        kl0 <- kl0 + xx.wt[i]*wk
    }
    list(ratio=kl/kl0,kl=kl)
}
