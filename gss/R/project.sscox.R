## Calculate Kullback-Leibler projection from sscox objects
project.sscox <- function(object,include,...)
{
    qd.pt <- object$mf
    if (is.null(object$cnt)) qd.wt <- rep(1,dim(qd.pt)[1])
    else qd.wt <- object$cnt
    if (!is.null(object$random))
        qd.wt <- qd.wt*exp(object$random$qd.z%*%object$b)
    bias <- object$bias
    ## evaluate full model
    mesh0 <- predict(object,qd.pt)
    qd.wt <- qd.wt*bias$qd.wt
    qd.wt <- t(t(qd.wt)/apply(qd.wt*mesh0,2,sum))
    ## extract terms in subspace
    nqd <- dim(qd.wt)[1]
    nxi <- length(object$id.basis)
    qd.s <- qd.r <- q <- NULL
    theta <- d <- NULL
    n0.wk <- nq.wk <- nq <- 0
    for (label in object$terms$labels) {
        x.basis <- object$mf[object$id.basis,object$term[[label]]$vlist]
        qd.x <- qd.pt[,object$term[[label]]$vlist]
        nphi <- object$term[[label]]$nphi
        nrk <- object$term[[label]]$nrk
        if (nphi) {
            phi <- object$term[[label]]$phi
            for (i in 1:nphi) {
                n0.wk <- n0.wk + 1
                if (!any(label==include)) next
                d <- c(d,object$d[n0.wk])
                qd.s <- cbind(qd.s,phi$fun(qd.x,nu=i,env=phi$env))
            }
        }
        if (nrk) {
            rk <- object$term[[label]]$rk
            for (i in 1:nrk) {
                nq.wk <- nq.wk + 1
                if (!any(label==include)) next
                nq <- nq + 1
                theta <- c(theta,object$theta[nq.wk])
                qd.r <- array(c(qd.r,rk$fun(x.basis,qd.x,nu=i,env=rk$env,out=TRUE)),
                              c(nxi,nqd,nq))
                q <- cbind(q,rk$fun(x.basis,x.basis,nu=i,env=rk$env,out=FALSE))
            }
        }
    }
    if (!is.null(object$partial)) {
        matx.p <- model.matrix(object$partial$mt,object$mf)[,-1,drop=FALSE]
        matx.p <- scale(matx.p)
        for (label in object$lab.p) {
            n0.wk <- n0.wk + 1
            if (!any(label==include)) next
            d <- c(d,object$d[n0.wk])
            qd.s <- cbind(qd.s,matx.p[,label])
        }
    }
    if (!is.null(qd.s)) {
        nn <- nxi + ncol(qd.s)
        qd.s <- t(qd.s)
    }
    else nn <- nxi
    ## calculate projection
    rkl <- function(theta1=NULL) {
        theta.wk <- 1:nq
        theta.wk[fix] <- theta[fix]
        if (nq-1) theta.wk[-fix] <- theta1
        qd.rs <- 0
        for (i in 1:nq) qd.rs <- qd.rs + 10^theta.wk[i]*qd.r[,,i]
        qd.rs <- rbind(qd.rs,qd.s)
        z <- .Fortran("drkl",
                      cd=as.double(cd), as.integer(nn),
                      as.double(t(qd.rs)), as.integer(nqd), as.integer(bias$nt),
                      as.double(bias$wt), as.double(t(qd.wt)),
                      mesh=as.double(mesh0), as.double(.Machine$double.eps),
                      as.double(1e-6), as.integer(30), double(nn),
                      double(2*bias$nt*(nqd+1)+nn*(2*nn+4)), info=integer(1),
                      PACKAGE="gss")
        if (z$info==1)
            stop("gss error in project.sscox: Newton iteration diverges")
        if (z$info==2)
            warning("gss warning in project.sscox: Newton iteration fails to converge")
        assign("cd",z$cd,inherits=TRUE)
        assign("mesh1",z$mesh,inherits=TRUE)
        sum(bias$wt*(apply(qd.wt*log(mesh0/mesh1)*mesh0,2,sum)+
                     log(apply(qd.wt*mesh1,2,sum))))
    }
    cv.wk <- function(theta) cv.scale*rkl(theta)+cv.shift
    if (nq) {
        ## initialization
        if (is.null(qd.s)) theta.wk <- 0
        else {
            qd.r.wk <- 0
            for (i in 1:nq) qd.r.wk <- qd.r.wk + 10^theta[i]*qd.r[,,i]
            vv.s <- vv.r <- 0
            for (i in 1:bias$nt) {
                mu.s <- apply(qd.wt[,i]*qd.s,2,sum)/sum(qd.wt[,i])
                v.s <- apply(qd.wt[,i]*qd.s^2,2,sum)/sum(qd.wt[,i])
                v.s <- v.s - mu.s^2
                mu.r <- apply(qd.wt[,i]*qd.r.wk,2,sum)/sum(qd.wt[,i])
                v.r <- apply(qd.wt[,i]*qd.r.wk^2,2,sum)/sum(qd.wt[,i])
                v.r <- v.r - mu.r^2
                vv.s <- vv.s + bias$wt[i]*v.s
                vv.r <- vv.r + bias$wt[i]*v.r
            }
            theta.wk <- log10(sum(vv.s)/(nn-nxi)/sum(vv.r)*nxi) / 2
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
        nn <- nrow(qd.s)
        z <- .Fortran("drkl",
                      cd=as.double(d), as.integer(nn),
                      as.double(qd.s), as.integer(nqd), as.integer(bias$nt),
                      as.double(bias$wt), as.double(t(qd.wt)),
                      mesh=as.double(mesh0), as.double(.Machine$double.eps),
                      as.double(1e-6), as.integer(30), double(nn),
                      double(2*bias$nt*(nqd+1)+nn*(2*nn+4)), info=integer(1),
                      PACKAGE="gss")
        if (z$info==1)
            stop("gss error in project.sscox: Newton iteration diverges")
        if (z$info==2)
            warning("gss warning in project.sscox: Newton iteration fails to converge")
        mesh1 <- z$mesh
        kl <- sum(bias$wt*(apply(qd.wt*log(mesh0/mesh1)*mesh0,2,sum)+
                           log(apply(qd.wt*mesh1,2,sum))))
    }
    kl0 <- sum(bias$wt*(apply(qd.wt*log(mesh0)*mesh0,2,sum)+
                        log(apply(qd.wt,2,sum))))
    wt.wk <- t(t(qd.wt)/apply(qd.wt*mesh1,2,sum))
    kl1 <- sum(bias$wt*(apply(wt.wk*log(mesh1)*mesh1,2,sum)+
                        log(apply(wt.wk,2,sum))))
    obj <- list(ratio=kl/kl0,kl=kl,check=(kl+kl1)/kl0)
    obj
}
