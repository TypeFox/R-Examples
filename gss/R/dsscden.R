dsscden <- ## Evaluate conditional density estimate
function (object,y,x) {
    ## check input
    if (!("sscden"%in%class(object))) stop("gss error in dsscden: not a sscden object")
    if (!all(sort(object$xnames)==sort(colnames(x))))
        stop("gss error in dsscden: mismatched x variable names")
    if (length(object$ynames)==1&is.vector(y)) {
        y <- data.frame(y)
        colnames(y) <- object$ynames
    }
    if (!all(sort(object$ynames)==sort(colnames(y))))
        stop("gss error in dsscden: mismatched y variable names")
    if ("sscden1"%in%class(object)) {
        qd.pt <- object$rho$env$qd.pt
        qd.wt <- object$rho$env$qd.wt
        d.qd <- d.sscden1(object,x,qd.pt,scale=FALSE)
        int <- apply(d.qd*qd.wt,2,sum)
        return(t(t(d.sscden1(object,x,y,scale=FALSE))/int))
    }
    else {
        qd.pt <- object$yquad$pt
        qd.wt <- object$yquad$wt
        d.qd <- d.sscden(object,x,qd.pt)
        int <- apply(d.qd*qd.wt,2,sum)
        return(t(t(d.sscden(object,x,y))/int))
    }
}

psscden <- ## Compute cdf for univariate density estimate
function(object,q,x) {
    if (!("sscden"%in%class(object))) stop("gss error in psscden: not a sscden object")
    if (length(object$ynames)!=1) stop("gss error in psscden: y is not 1-D")
    if (("sscden1"%in%class(object))&!is.numeric(object$mf[,object$ynames]))
        stop("gss error in qssden: y is not continuous")
    if ("sscden1"%in%class(object)) ydomain <- object$rho$env$ydomain
    else ydomain <- object$ydomain
    mn <- min(ydomain[[object$ynames]])
    mx <- max(ydomain[[object$ynames]])
    order.q <- rank(q)
    p <- q <- sort(q)
    q.dup <- duplicated(q)
    p[q<=mn] <- 0
    p[q>=mx] <- 1
    qd.hize <- 200
    qd <- gauss.quad(2*qd.hize,c(mn,mx))
    y.wk <- data.frame(qd$pt)
    colnames(y.wk) <- object$ynames
    d.qd <- dsscden(object,y.wk,x)
    gap <- diff(qd$pt)
    g.wk <- gap[qd.hize]/2
    for (i in 1:(qd.hize-2)) g.wk <- c(g.wk,gap[qd.hize+i]-g.wk[i])
    g.wk <- 2*g.wk
    g.wk <- c(g.wk,(mx-mn)/2-sum(g.wk))
    gap[qd.hize:1] <- gap[qd.hize+(1:qd.hize)] <- g.wk
    brk <- cumsum(c(mn,gap))
    kk <- (1:length(q))[q>mn&q<mx]
    z <- NULL
    for (k in 1:dim(x)[1]) {
        d.qd.wk <- d.qd[,k]/sum(d.qd[,k]*qd$wt)
        for (i in kk) {
            if (q.dup[i]) {
                p[i] <- p.dup
                next
            }
            ind <- (1:(2*qd.hize))[qd$pt<q[i]]
            if (!length(ind)) {
                wk <- d.qd.wk[1]*qd$wt[1]*(q[i]-mn)/gap[1]
            }
            else {
                wk <- sum(d.qd.wk[ind]*qd$wt[ind])
                id.mx <- max(ind)
                if (q[i]<brk[id.mx+1])
                  wk <- wk-d.qd.wk[id.mx]*qd$wt[id.mx]*(brk[id.mx+1]-q[i])/gap[id.mx]
                else wk <- wk+d.qd.wk[id.mx+1]*qd$wt[id.mx+1]*(q[i]-brk[id.mx+1])/gap[id.mx+1]
            }
            p[i] <- p.dup <- wk
        }
        z <- cbind(z,p[order.q])
    }
    z
}

qsscden <- ## Compute quantiles for univariate density estimate
function(object,p,x) {
    if (!("sscden"%in%class(object))) stop("gss error in qsscden: not a sscden object")
    if (length(object$ynames)!=1) stop("gss error in qsscden: y is not 1-D")
    if (("sscden1"%in%class(object))&!is.numeric(object$mf[,object$ynames]))
        stop("gss error in qssden: y is not continuous")
    if ("sscden1"%in%class(object)) ydomain <- object$rho$env$ydomain
    else ydomain <- object$ydomain
    mn <- min(ydomain[[object$ynames]])
    mx <- max(ydomain[[object$ynames]])
    order.p <- rank(p)
    q <- p <- sort(p)
    p.dup <- duplicated(p)
    q[p<=0] <- mn
    q[p>=1] <- mx
    qd.hize <- 200
    qd <- gauss.quad(2*qd.hize,c(mn,mx))
    y.wk <- data.frame(qd$pt)
    colnames(y.wk) <- object$ynames
    d.qd <- dsscden(object,y.wk,x)
    gap <- diff(qd$pt)
    g.wk <- gap[qd.hize]/2
    for (i in 1:(qd.hize-2)) g.wk <- c(g.wk,gap[qd.hize+i]-g.wk[i])
    g.wk <- 2*g.wk
    g.wk <- c(g.wk,(mx-mn)/2-sum(g.wk))
    gap[qd.hize:1] <- gap[qd.hize+(1:qd.hize)] <- g.wk
    brk <- cumsum(c(mn,gap))
    kk <- (1:length(p))[p>0&p<1]
    z <- NULL
    for (k in 1:dim(x)[1]) {
        d.qd.wk <- d.qd[,k]/sum(d.qd[,k]*qd$wt)
        p.wk <- cumsum(d.qd.wk*qd$wt)
        for (i in kk) {
            if (p.dup[i]) {
                q[i] <- q.dup
                next
            }
            ind <- (1:(2*qd.hize))[p.wk<p[i]]
            if (!length(ind)) {
                wk <- mn+p[i]/(d.qd.wk[1]*qd$wt[1])*gap[1]
            }
            else {
                id.mx <- max(ind)
                wk <- brk[id.mx+1]+(p[i]-p.wk[id.mx])/(d.qd.wk[id.mx+1]*qd$wt[id.mx+1])*gap[id.mx+1]
            }
            q[i] <- q.dup <- wk
        }
        z <- cbind(z,q[order.p])
    }
    z
}

d.sscden <- ## Evaluate conditional density estimate
function (object,x,y) {
    ## check input
    if ("sscden"!=class(object)) stop("gss error in d.sscden: not a sscden object")
    if (!all(sort(object$xnames)==sort(colnames(x))))
        stop("gss error in d.sscden: mismatched x variable names")
    if (!all(sort(object$ynames)==sort(colnames(y))))
        stop("gss error in d.sscden: mismatched y variable names")
    mf <- object$mf
    ## exp(eta)
    z <- NULL
    for (k in 1:dim(x)[1]) {
        xy.new <- cbind(x[rep(k,dim(y)[1]),,drop=FALSE],y)
        s <- NULL
        r <- 0
        nu <- nq <- 0
        for (label in object$terms$labels) {
            vlist <- object$terms[[label]]$vlist
            xy <- xy.new[,vlist]
            xy.basis <- mf[object$id.basis,vlist]
            nphi <- object$terms[[label]]$nphi
            nrk <- object$terms[[label]]$nrk
            if (nphi) {
                phi <- object$terms[[label]]$phi
                for (i in 1:nphi) {
                    nu <- nu + 1
                    s.wk <- phi$fun(xy,nu=i,env=phi$env)
                    s <- cbind(s,s.wk)
                }
            }
            if (nrk) {
                rk <- object$terms[[label]]$rk
                for (i in 1:nrk) {
                    nq <- nq+1
                    r.wk <- rk$fun(xy,xy.basis,nu=i,env=rk$env,out=TRUE)
                    r <- r + 10^object$theta[nq]*r.wk
                }
            }
        }
        eta <- exp(cbind(s,r)%*%c(object$d,object$c))
        z <- cbind(z,eta)
    }
    z
}

d.sscden1 <- ## Evaluate conditional density estimate
function (object,x,y,scale=TRUE) {
    ## check input
    if (!("sscden1"%in%class(object))) stop("gss error in d.sscden1: not a sscden1 object")
    if (!all(sort(object$xnames)==sort(colnames(x))))
        stop("gss error in d.sscden1: mismatched x variable names")
    if (!all(sort(object$ynames)==sort(colnames(y))))
        stop("gss error in d.sscden1: mismatched y variable names")
    mf <- object$mf
    ## rho
    rho <- object$rho$fun(x,y,object$rho$env,outer=TRUE)
    ## exp(eta)
    z <- NULL
    for (k in 1:dim(x)[1]) {
        xy.new <- cbind(x[rep(k,dim(y)[1]),,drop=FALSE],y)
        s <- NULL
        r <- 0
        nu <- nq <- 0
        for (label in object$terms$labels) {
            vlist <- object$terms[[label]]$vlist
            xy <- xy.new[,vlist]
            xy.basis <- mf[object$id.basis,vlist]
            nphi <- object$terms[[label]]$nphi
            nrk <- object$terms[[label]]$nrk
            if (nphi) {
                phi <- object$terms[[label]]$phi
                for (i in 1:nphi) {
                    nu <- nu + 1
                    if (!scale&!(nu%in%object$id.s)) next
                    s.wk <- phi$fun(xy,nu=i,env=phi$env)
                    s <- cbind(s,s.wk)
                }
            }
            if (nrk) {
                rk <- object$terms[[label]]$rk
                for (i in 1:nrk) {
                    nq <- nq+1
                    if (!scale&!(nq%in%object$id.r)) next
                    r.wk <- rk$fun(xy,xy.basis,nu=i,env=rk$env,out=TRUE)
                    r <- r + 10^object$theta[nq]*r.wk
                }
            }
        }
        if (!scale) eta <- exp(cbind(s,r)%*%c(object$d[object$id.s],object$c))
        else eta <- exp(cbind(s,r)%*%c(object$d,object$c))*object$scal
        z <- cbind(z,eta*rho[k,])
    }
    z
}
