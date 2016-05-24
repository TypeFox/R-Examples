## Calculate square error projection from sscden objects
project.sscden1 <- function(object,include,...)
{
    include <- unique(include)
    term <- object$terms
    mf <- object$mf
    nobs <- dim(mf)[1]
    xnames <- object$xnames
    ynames <- object$ynames
    xx.wt <- object$xx.wt
    nx <- length(xx.wt)
    xx <- mf[!object$x.dup.ind,xnames,drop=FALSE]
    qd.pt <- object$rho$env$qd.pt
    qd.wt <- object$rho$env$qd.wt
    rho.d <- t(object$rho$fun(xx,qd.pt,object$rho$env,outer=TRUE))
    rho.wt <- rho.d*qd.wt
    rho.d <- t(t(log(rho.d))-apply(log(rho.d)*rho.wt,2,sum))
    nmesh <- length(qd.wt)
    ns <- length(object$id.s)
    nr <- length(object$id.r)
    nbasis <- length(object$id.basis)
    s.rho <- ss <- 0
    r.rho <- matrix(0,nbasis,nr)
    sr <- array(0,c(ns,nbasis,nr))
    rr <- array(0,c(nbasis,nbasis,nr,nr))
    rho2 <- sum(xx.wt*apply(rho.d^2*rho.wt,2,sum))
    for (k in 1:nx) {
        id.x <- (1:nobs)[!object$x.dup.ind][k]
        qd.s <- NULL
        qd.r <- list(NULL)
        iq <- 0
        for (label in term$labels) {
            vlist <- term[[label]]$vlist
            x.list <- xnames[xnames%in%vlist]
            y.list <- ynames[ynames%in%vlist]
            xy.basis <- mf[object$id.basis,vlist]
            qd.xy <- data.frame(matrix(0,nmesh,length(vlist)))
            names(qd.xy) <- vlist
            qd.xy[,y.list] <- qd.pt[,y.list]
            if (length(x.list)) xx <- mf[rep(id.x,nmesh),x.list,drop=FALSE]
            else xx <- NULL
            nphi <- term[[label]]$nphi
            nrk <- term[[label]]$nrk
            if (nphi) {
                phi <- term[[label]]$phi
                for (i in 1:nphi) {
                    if (is.null(xx)) {
                        qd.s.wk <- phi$fun(qd.xy[,,drop=TRUE],nu=i,env=phi$env)
                        qd.s <- cbind(qd.s,qd.s.wk)
                    }
                    else {
                        if (length(y.list)>0) {
                            qd.xy[,x.list] <- xx
                            qd.s <- cbind(qd.s,phi$fun(qd.xy,i,phi$env))
                        }
                    }
                }
            }
            if (nrk) {
                rk <- term[[label]]$rk
                for (i in 1:nrk) {
                    if (is.null(xx)) {
                        qd.r.wk <- rk$fun(qd.xy[,,drop=TRUE],xy.basis,nu=i,env=rk$env,TRUE)
                        iq <- iq+1
                        qd.r[[iq]] <- qd.r.wk
                    }
                    else {
                        if (length(y.list)>0) {
                            qd.xy[,x.list] <- xx
                            iq <- iq+1
                            qd.r[[iq]] <- rk$fun(qd.xy,xy.basis,i,rk$env,TRUE)
                        }
                    }
                }
            }
        }
        if (ns) {
            qd.s <- sweep(qd.s,2,apply(qd.s*rho.wt[,k],2,sum))
            s.rho <- s.rho + xx.wt[k]*apply(qd.s*rho.d[,k]*rho.wt[,k],2,sum)
            ss <- ss + xx.wt[k]*t(rho.wt[,k]*qd.s)%*%qd.s
        }
        for (i in 1:iq) {
            qd.r[[i]] <- sweep(qd.r[[i]],2,apply(qd.r[[i]]*rho.wt[,k],2,sum))
            r.rho[,i] <- r.rho[,i] + xx.wt[k]*apply(qd.r[[i]]*rho.d[,k]*rho.wt[,k],2,sum)
            if (ns) sr[,,i] <- sr[,,i] + xx.wt[k]*t(rho.wt[,k]*qd.s)%*%qd.r[[i]]
            for (j in 1:i) {
                rr.wk <- xx.wt[k]*t(rho.wt[,k]*qd.r[[i]])%*%qd.r[[j]]
                rr[,,i,j] <- rr[,,i,j] + rr.wk
                if (i-j) rr[,,j,i] <- rr[,,j,i] + t(rr.wk)
            }
        }
    }
    ## evaluate full model
    if (ns) d <- object$d[object$id.s]
    c <- object$c
    theta <- object$theta[object$id.r]
    nq <- length(theta)
    if (ns) s.eta <- ss%*%d
    r.eta <- tmp <- NULL
    r.rho.wk <- sr.wk <- rr.wk <- 0
    for (i in 1:nq) {
        tmp <- c(tmp,10^(2*theta[i])*sum(diag(rr[,,i,i])))
        if (ns) {
            s.eta <- s.eta + 10^theta[i]*sr[,,i]%*%c
            if (length(d)==1) r.eta.wk <- sr[,,i]*d
            else r.eta.wk <- t(sr[,,i])%*%d
            sr.wk <- sr.wk + 10^theta[i]*sr[,,i]
        }
        else r.eta.wk <- 0
        r.rho.wk <- r.rho.wk + 10^theta[i]*r.rho[,i]
        for (j in 1:nq) {
            r.eta.wk <- r.eta.wk + 10^theta[j]*rr[,,i,j]%*%c
            rr.wk <- rr.wk + 10^(theta[i]+theta[j])*rr[,,i,j]
        }
        r.eta <- cbind(r.eta,r.eta.wk)
    }
    rho.eta <- sum(r.rho.wk*c)
    if (ns) rho.eta <- rho.eta + sum(r.rho.wk*c)
    eta2 <- sum(c*(rr.wk%*%c))
    if (ns) eta2 <- eta2 + sum(d*(ss%*%d)) + 2*sum(d*(sr.wk%*%c))
    mse <- eta2 + rho2 + 2*rho.eta
    ## extract terms in subspace
    id.s <- id.r <- NULL
    for (label in include) {
        id.s <- c(id.s,object$id.s.list[[label]])
        id.r <- c(id.r,object$id.r.list[[label]])
    }
    if (is.null(id.s)&is.null(id.r))
        stop("gss error in project.sscden1: include some terms")
    if (!all(id.s%in%object$id.s)|!all(id.r%in%object$id.r))
        stop("gss error in project.sscden1: included terms are not in the model")
    ## calculate projection
    rkl <- function(theta1=NULL) {
        theta.wk <- 1:nq0
        theta.wk[fix] <- theta[fix]
        if (nq0-1) theta.wk[-fix] <- theta1
        ##
        id.s0 <- (1:length(object$id.s))[object$id.s%in%id.s]
        id.r0 <- (1:length(object$id.r))[object$id.r%in%id.r]
        if (length(id.s0)) ss.wk <- ss[id.s0,id.s0,drop=FALSE]
        if (length(id.r0)) {
            r.eta.wk <- rr.wk <- 0
            sr.wk <- matrix(0,length(id.s),nbasis)
            for (i in 1:length(id.r0)) {
                r.eta.wk <- r.eta.wk + 10^theta.wk[i]*r.eta[,id.r0[i]]
                if (length(id.s0)) sr.wk <- sr.wk + 10^theta.wk[i]*sr[id.s0,,id.r0[i]]
                for (j in 1:length(id.r0)) {
                    rr.wk <- rr.wk + 10^(theta.wk[i]+theta.wk[j])*rr[,,id.r0[i],id.r0[j]]
                }
            }
            if (length(id.s0)) {
                v <- cbind(rbind(ss.wk,t(sr.wk)),rbind(sr.wk,rr.wk))
                mu <- c(s.eta[id.s0],r.eta.wk)
            }
            else {
                v <- rbind(sr.wk,rr.wk)
                mu <- r.eta.wk
            }
        }
        else {
            v <- ss.wk
            mu <- s.eta[id.s0]
        }
        nn <- length(mu)
        z <- chol(v,pivot=TRUE)
        v <- z
        rkv <- attr(z,"rank")
        m.eps <- .Machine$double.eps
        while (v[rkv,rkv]<2*sqrt(m.eps)*v[1,1]) rkv <- rkv - 1
        if (rkv<nn) v[(1:nn)>rkv,(1:nn)>rkv] <- diag(v[1,1],nn-rkv)
        mu <- backsolve(v,mu[attr(z,"pivot")],transpose=TRUE)
        eta2 - sum(mu[1:rkv]^2)
    }
    cv.wk <- function(theta) cv.scale*rkl(theta)+cv.shift
    ## initialization
    fix <- rev(order(tmp[id.r]))[1]
    theta <- object$theta[id.r]
    ## projection
    nq0 <- length(id.r)
    if (nq0>1) {
        if (object$skip.iter) se <- rkl(theta[-fix])
        else {
            if (nq0-2) {
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
            se <- rkl(zz$est)
        }
    }
    else se <- rkl()
    list(ratio=se/mse,se=se)
}
