## Calculate square error projection from sshzd1 objects
project.sshzd1 <- function(object,include,...)
{
    if (!(object$tname%in%include))
        stop("gss error in project.sshzd1: time main effect missing in included terms")
    ## Initialization
    term <- object$term
    mf <- object$mf
    xnames <- object$xnames
    tname <- object$tname
    id.basis <- object$id.basis
    yy <- object$yy
    quad <- object$quad
    x.pt <- object$x.pt
    qd.wt <- object$qd.wt
    ## Calculate cross integrals of phi and rk
    s <- object$int.s
    r <- object$int.r
    ns <- length(s)
    nq <- length(object$theta)
    nx <- dim(qd.wt)[2]
    nbasis <- dim(r)[1]
    ## create arrays
    ss <- 0
    sr <- array(0,c(ns,nbasis,nq))
    rr <- array(0,c(nbasis,nbasis,nq,nq))
    for (k in 1:nx) {
        ind <- (1:length(quad$pt))[qd.wt[,k]>0]
        nmesh <- length(ind)
        if (!nmesh) next
        qd.wt.wk <- qd.wt[ind,k]
        qd.s <- NULL
        qd.r <- as.list(NULL)
        iq <- 0
        for (label in term$labels) {
            if (label=="1") {
                qd.wk <- rep(1,nmesh)
                qd.s <- cbind(qd.s,qd.wk)
                next
            }
            vlist <- term[[label]]$vlist
            x.list <- xnames[xnames%in%vlist]
            xy.basis <- mf[id.basis,vlist]
            qd.xy <- data.frame(matrix(0,nmesh,length(vlist)))
            names(qd.xy) <- vlist
            if (tname%in%vlist) qd.xy[,tname] <- quad$pt[ind]
            if (length(x.list)) qd.xy[,x.list] <- x.pt[rep(k,nmesh),x.list,drop=FALSE]
            nphi <- term[[label]]$nphi
            nrk <- term[[label]]$nrk
            if (nphi) {
                phi <- term[[label]]$phi
                for (i in 1:nphi) {
                    qd.wk <- phi$fun(qd.xy[,,drop=TRUE],nu=i,env=phi$env)
                    qd.s <- cbind(qd.s,qd.wk)
                }
            }
            if (nrk) {
                rk <- term[[label]]$rk
                for (i in 1:nrk) {
                    iq <- iq+1
                    qd.r[[iq]] <- rk$fun(qd.xy[,,drop=TRUE],xy.basis,i,rk$env,out=TRUE)
                }
            }
        }
        if (!is.null(object$partial)) {
            wk <- object$partial$pt[k,]
            qd.s <- cbind(qd.s,t(matrix(wk,length(wk),nmesh)))
        }
        ss <- ss + t(qd.wt.wk*qd.s)%*%qd.s
        for (i in 1:nq) {
            sr[,,i] <- sr[,,i] + t(qd.wt.wk*qd.s)%*%qd.r[[i]]
            for (j in 1:i) {
                rr.wk <- t(qd.wt.wk*qd.r[[i]])%*%qd.r[[j]]
                rr[,,i,j] <- rr[,,i,j] + rr.wk
                if (i-j) rr[,,j,i] <- rr[,,j,i] + t(rr.wk)
            }
        }
    }
    ## evaluate full model
    cfit <- log(object$cfit)
    d <- object$d
    c <- object$c
    theta <- object$theta
    s.eta <- ss%*%d
    r.eta <- tmp <- NULL
    r.wk <- sr.wk <- rr.wk <- 0
    for (i in 1:nq) {
        tmp <- c(tmp,10^(2*theta[i])*sum(diag(rr[,,i,i])))
        s.eta <- s.eta + 10^theta[i]*sr[,,i]%*%c
        if (length(d)==1) r.eta.wk <- sr[,,i]*d
        else r.eta.wk <- t(sr[,,i])%*%d
        r.wk <- r.wk + 10^theta[i]*r[,i]
        sr.wk <- sr.wk + 10^theta[i]*sr[,,i]
        for (j in 1:nq) {
            r.eta.wk <- r.eta.wk + 10^theta[j]*rr[,,i,j]%*%c
            rr.wk <- rr.wk + 10^(theta[i]+theta[j])*rr[,,i,j]
        }
        r.eta <- cbind(r.eta,r.eta.wk)
    }
    eta2 <- sum(c*(rr.wk%*%c)) + sum(d*(ss%*%d)) + 2*sum(d*(sr.wk%*%c))
    mse <- eta2 - 2*sum(c(d,c)*c(s,r.wk))*cfit + cfit^2*sum(qd.wt)
    ## extract terms in subspace
    id.s <- id.q <- NULL
    for (label in term$labels) {
        if (label=="1") {
            id.s <- c(id.s,1)
            next
        }
        if (!any(label==include)) next
        term.wk <- term[[label]]
        if (term.wk$nphi>0) id.s <- c(id.s,term.wk$iphi+(1:term.wk$nphi)-1)
        if (term.wk$nrk>0) id.q <- c(id.q,term.wk$irk+(1:term.wk$nrk)-1)
    }
    if (!is.null(object$partial)) {
        nu <- length(object$d)-length(object$lab.p)
        for (label in object$lab.p) {
            nu <- nu+1
            if (!any(label==include)) next
            id.s <- c(id.s,nu)
        }
    }
    ## calculate projection
    rkl <- function(theta1=NULL) {
        theta.wk <- 1:nq
        theta.wk[fix] <- theta[fix]
        if (nq0-1) theta.wk[id.q0] <- theta1
        ##
        ss.wk <- ss[id.s,id.s]
        r.eta.wk <- sr.wk <- rr.wk <- 0
        for (i in id.q) {
            r.eta.wk <- r.eta.wk + 10^theta.wk[i]*r.eta[,i]
            sr.wk <- sr.wk + 10^theta.wk[i]*sr[id.s,,i]
            for (j in id.q) {
                rr.wk <- rr.wk + 10^(theta.wk[i]+theta.wk[j])*rr[,,i,j]
            }
        }
        v <- cbind(rbind(ss.wk,t(sr.wk)),rbind(sr.wk,rr.wk))
        mu <- c(s.eta[id.s],r.eta.wk)
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
    nq0 <- length(id.q)
    tmp[-id.q] <- 0
    fix <- rev(order(tmp))[1]
    ## projection
    if (nq0-1) {
        id.q0 <- id.q[id.q!=fix]
        if (object$skip.iter) se <- rkl(theta[id.q0])
        else {
            if (nq0-2) {
                ## scale and shift cv
                tmp <- abs(rkl(theta[id.q0]))
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
                zz <- nlm(cv.wk,theta[id.q0],stepmax=.5,ndigit=7)
            }
            else {
                the.wk <- theta[id.q0]
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
