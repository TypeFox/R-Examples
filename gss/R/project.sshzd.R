## Calculate Kullback-Leibler projection from sshzd objects
project.sshzd <- function(object,include,mesh=FALSE,...)
{
    if (!(object$tname%in%include))
        stop("gss error in project.sshzd: time main effect missing in included terms")
    quad.pt <- object$quad$pt
    qd.wt <- object$qd.wt
    nx <- dim(object$qd.wt)[2]
    nbasis <- length(object$id.basis)
    mesh0 <- object$mesh0
    ## extract terms in subspace
    nqd <- length(quad.pt)
    nxi <- length(object$id.basis)
    d <- qd.s <- q <- theta <- NULL
    qd.r <- as.list(NULL)
    n0.wk <- nu <- nq.wk <- nq <- 0
    for (label in object$terms$labels) {
        vlist <- object$terms[[label]]$vlist
        x.list <- object$xnames[object$xnames%in%vlist]
        xy.basis <- object$mf[object$id.basis,vlist]
        qd.xy <- data.frame(matrix(0,nqd,length(vlist)))
        names(qd.xy) <- vlist
        if (object$tname%in%vlist) qd.xy[,object$tname] <- quad.pt
        if (length(x.list)) xx <- object$x.pt[,x.list,drop=FALSE]
        else xx <- NULL
        nphi <- object$terms[[label]]$nphi
        nrk <- object$terms[[label]]$nrk
        if (nphi) {
            phi <- object$terms[[label]]$phi
            for (i in 1:nphi) {
                n0.wk <- n0.wk + 1
                if (label=="1") {
                    d <- object$d[n0.wk]
                    nu <- nu + 1
                    qd.wk <- matrix(1,nqd,nx)
                    qd.s <- array(c(qd.s,qd.wk),c(nqd,nx,nu))
                    next
                }
                if (!any(label==include)) next
                d <- c(d,object$d[n0.wk])
                nu <- nu + 1
                if (is.null(xx))
                    qd.wk <- matrix(phi$fun(qd.xy[,,drop=TRUE],nu=i,env=phi$env),nqd,nx)
                else {
                    qd.wk <- NULL
                    for (j in 1:nx) {
                        qd.xy[,x.list] <- xx[rep(j,nqd),]
                        for (k in x.list)
                            if (is.factor(xx[,k])) qd.xy[,k] <- as.factor(qd.xy[,k])
                        qd.wk <- cbind(qd.wk,phi$fun(qd.xy[,,drop=TRUE],i,phi$env))
                    }
                }
                qd.s <- array(c(qd.s,qd.wk),c(nqd,nx,nu))
            }
        }
        if (nrk) {
            rk <- object$terms[[label]]$rk
            for (i in 1:nrk) {
                nq.wk <- nq.wk + 1
                if (!any(label==include)) next
                nq <- nq + 1
                theta <- c(theta,object$theta[nq.wk])
                q <- cbind(q,rk$fun(xy.basis,xy.basis,i,rk$env,out=FALSE))
                if (is.null(xx))
                    qd.r[[nq]] <- rk$fun(qd.xy[,,drop=TRUE],xy.basis,i,rk$env,out=TRUE)
                else {
                    qd.wk <- NULL
                    for (j in 1:nx) {
                        qd.xy[,x.list] <- xx[rep(j,nqd),]
                        for (k in x.list)
                            if (is.factor(xx[,k])) qd.xy[,k] <- as.factor(qd.xy[,k])
                        qd.wk <- array(c(qd.wk,rk$fun(qd.xy[,,drop=TRUE],xy.basis,i,rk$env,TRUE)),
                                       c(nqd,nbasis,j))
                    }
                    qd.r[[nq]] <- qd.wk
                }
            }
        }
    }
    if (!is.null(object$partial)) {
        for (label in object$lab.p) {
            n0.wk <- n0.wk + 1
            if (!any(label==include)) next
            d <- c(d,object$d[n0.wk])
            qd.wk <- t(matrix(object$partial$pt[,label],nx,nqd))
            qd.s <- array(c(qd.s,qd.wk),c(nqd,nx,n0.wk))
        }
    }
    if (!is.null(qd.s)) nnull <- dim(qd.s)[3]
    else nnull <- 0
    nn <- nxi + nnull
    ## random effect offset
    if (!is.null(object$b)) offset <- as.vector(object$random$qd.z%*%object$b)
    else offset <- rep(0,nx)
    ## calculate projection
    rkl <- function(theta1=NULL) {
        theta.wk <- 1:nq
        theta.wk[fix] <- theta[fix]
        if (nq-1) theta.wk[-fix] <- theta1
        qd.r.wk <- array(0,c(nqd,nxi,nx))
        for (i in 1:nq) {
            if (length(dim(qd.r[[i]]))==3) qd.r.wk <- qd.r.wk + 10^theta.wk[i]*qd.r[[i]]
            else qd.r.wk <- qd.r.wk + as.vector(10^theta.wk[i]*qd.r[[i]])
        }
        qd.r.wk <- aperm(qd.r.wk,c(1,3,2))
        qd.r.wk <- array(c(qd.r.wk,qd.s),c(nqd,nx,nn))
        qd.r.wk <- aperm(qd.r.wk,c(1,3,2))
        z <- .Fortran("hrkl",
                      cd=as.double(cd), as.integer(nn),
                      as.double(qd.r.wk), as.integer(nqd), as.integer(nx),
                      as.double(t(t(qd.wt)*exp(offset))),
                      mesh=as.double(qd.wt*mesh0),
                      as.double(.Machine$double.eps), double(nqd*nx),
                      double(nn), double(nn), double(nn*nn), integer(nn), double(nn),
                      double(nn), double(nqd*nx), as.double(1e-6), as.integer(30),
                      info=integer(1), PACKAGE="gss")
        if (z$info==1)
            stop("gss error in project.sshzd: Newton iteration diverges")
        if (z$info==2)
            warning("gss warning in project.sshzd: Newton iteration fails to converge")
        assign("cd",z$cd,inherits=TRUE)
        assign("mesh1",t(t(matrix(z$mesh,nqd,nx))*exp(offset)),inherits=TRUE)
        sum(qd.wt*(log(mesh0/mesh1)*mesh0-mesh0+mesh1))
    }
    cv.wk <- function(theta) cv.scale*rkl(theta)+cv.shift
    ## initialization
    if (nnull) {
        qd.r.wk <- array(0,c(nqd,nxi,nx))
        for (i in 1:nq) {
            if (length(dim(qd.r[[i]]))==3) qd.r.wk <- qd.r.wk + 10^theta[i]*qd.r[[i]]
            else qd.r.wk <- qd.r.wk + as.vector(10^theta[i]*qd.r[[i]])
        }
        v.s <- v.r <- 0
        for (i in 1:nx) {
            v.s <- v.s + apply(qd.wt[,i]*qd.s[,i,,drop=FALSE]^2,2,sum)
            v.r <- v.r + apply(qd.wt[,i]*qd.r.wk[,,i,drop=FALSE]^2,2,sum)
        }
        theta.wk <- log10(sum(v.s)/nnull/sum(v.r)*nxi) / 2
    }
    else theta.wk <- 0
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
    ## cfit
    cfit <- t(matrix(object$dbar/sum(t(qd.wt)*exp(offset))*exp(offset),nx,nqd))
    ## return
    kl0 <- sum(object$qd.wt*(log(mesh0/cfit)*mesh0-mesh0+cfit))
    kl1 <- sum(object$qd.wt*(log(mesh1/cfit)*mesh1-mesh1+cfit))
    obj <- list(ratio=kl/kl0,kl=kl,check=(kl+kl1)/kl0)
    if (mesh) obj$mesh <- mesh1
    obj
}
