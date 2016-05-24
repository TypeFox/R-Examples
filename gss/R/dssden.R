dssden <- ## Evaluate density estimate
function (object,x) {
    ## check input
    if (!("ssden"%in%class(object))) stop("gss error in dssden: not a ssden object")
    if ("ssden1"%in%class(object)) return(d.ssden1(object,x))
    else return(d.ssden(object,x))
}

pssden <- ## Compute cdf for univariate density estimate
function(object,q) {
    if (!("ssden"%in%class(object))) stop("gss error in pssden: not a ssden object")
    if (dim(object$mf)[2]!=1) stop("gss error in pssden: not a 1-D density")
    mn <- min(object$domain)
    mx <- max(object$domain)
    order.q <- rank(q)
    p <- q <- sort(q)
    q.dup <- duplicated(q)
    p[q<=mn] <- 0
    p[q>=mx] <- 1
    qd.hize <- 200
    qd <- gauss.quad(2*qd.hize,c(mn,mx))
    d.qd <- dssden(object,qd$pt)
    gap <- diff(qd$pt)
    g.wk <- gap[qd.hize]/2
    for (i in 1:(qd.hize-2)) g.wk <- c(g.wk,gap[qd.hize+i]-g.wk[i])
    g.wk <- 2*g.wk
    g.wk <- c(g.wk,(mx-mn)/2-sum(g.wk))
    gap[qd.hize:1] <- gap[qd.hize+(1:qd.hize)] <- g.wk
    brk <- cumsum(c(mn,gap))
    kk <- (1:length(q))[q>mn&q<mx]
    for (i in kk) {
        if (q.dup[i]) {
            p[i] <- p.dup
            next
        }
        ind <- (1:(2*qd.hize))[qd$pt<q[i]]
        if (!length(ind)) {
            wk <- d.qd[1]*qd$wt[1]*(q[i]-mn)/gap[1]
        }
        else {
            wk <- sum(d.qd[ind]*qd$wt[ind])
            id.mx <- max(ind)
            if (q[i]<brk[id.mx+1])
                wk <- wk-d.qd[id.mx]*qd$wt[id.mx]*(brk[id.mx+1]-q[i])/gap[id.mx]
            else wk <- wk+d.qd[id.mx+1]*qd$wt[id.mx+1]*(q[i]-brk[id.mx+1])/gap[id.mx+1]
        }
        p[i] <- p.dup <- wk
    }
    p[order.q]
}

qssden <- ## Compute quantiles for univariate density estimate
function(object,p) {
    if (!("ssden"%in%class(object))) stop("gss error in qssden: not a ssden object")
    if (dim(object$mf)[2]!=1) stop("gss error in qssden: not a 1-D density")
    mn <- min(object$domain)
    mx <- max(object$domain)
    order.p <- rank(p)
    q <- p <- sort(p)
    p.dup <- duplicated(p)
    q[p<=0] <- mn
    q[p>=1] <- mx
    qd.hize <- 200
    qd <- gauss.quad(2*qd.hize,c(mn,mx))
    d.qd <- dssden(object,qd$pt)
    gap <- diff(qd$pt)
    g.wk <- gap[qd.hize]/2
    for (i in 1:(qd.hize-2)) g.wk <- c(g.wk,gap[qd.hize+i]-g.wk[i])
    g.wk <- 2*g.wk
    g.wk <- c(g.wk,(mx-mn)/2-sum(g.wk))
    gap[qd.hize:1] <- gap[qd.hize+(1:qd.hize)] <- g.wk
    brk <- cumsum(c(mn,gap))
    p.wk <- cumsum(d.qd*qd$wt)
    kk <- (1:length(p))[p>0&p<1]
    for (i in kk) {
        if (p.dup[i]) {
            q[i] <- q.dup
            next
        }
        ind <- (1:(2*qd.hize))[p.wk<p[i]]
        if (!length(ind)) {
            wk <- mn+p[i]/(d.qd[1]*qd$wt[1])*gap[1]
        }
        else {
            id.mx <- max(ind)
            wk <- brk[id.mx+1]+(p[i]-p.wk[id.mx])/(d.qd[id.mx+1]*qd$wt[id.mx+1])*gap[id.mx+1]
        }
        q[i] <- q.dup <- wk
    }
    q[order.p]
}

d.ssden <- ## Evaluate density estimate
function (object,x) {
    if (class(object)!="ssden") stop("gss error in d.ssden: not a ssden object")
    if (dim(object$mf)[2]==1&is.vector(x)) {
        x <- data.frame(x)
        colnames(x) <- colnames(object$mf)
    }
    s <- NULL
    r <- matrix(0,dim(x)[1],length(object$id.basis))
    nq <- 0
    for (label in object$terms$labels) {
        xx <- object$mf[object$id.basis,object$terms[[label]]$vlist]
        x.new <- x[,object$terms[[label]]$vlist]
        nphi <- object$terms[[label]]$nphi
        nrk <-  object$terms[[label]]$nrk
        if (nphi) {
            phi <-  object$terms[[label]]$phi
            for (i in 1:nphi) {
                s <- cbind(s,phi$fun(x.new,nu=i,env=phi$env))
            }
        }
        if (nrk) {
            rk <- object$terms[[label]]$rk
            for (i in 1:nrk) {
                nq <- nq + 1
                r <- r + 10^object$theta[nq]*rk$fun(x.new,xx,nu=i,env=rk$env,out=TRUE)
            }
        }
    }
    as.vector(exp(cbind(s,r)%*%c(object$d,object$c))/object$int)
}

d.ssden1 <- ## Evaluate density estimate
function (object,x) {
    if (!("ssden1"%in%class(object))) stop("gss error in d.ssden1: not a ssden1 object")
    ## rho
    rho <- 1
    for (xlab in names(object$mf)) {
        xx <- x[[xlab]]
        rho.wk <- object$rho[[xlab]]
        if (is.factor(xx)) rho <- rho*rho.wk[xx]
        if (is.vector(xx)&!is.factor(xx)) rho <- rho*dssden(rho.wk,xx)
        if (is.matrix(xx)) rho <- rho*dssden(rho.wk,xx)
    }
    ## exp(eta)
    s <- NULL
    r <- matrix(0,dim(x)[1],length(object$id.basis))
    nq <- 0
    for (label in object$terms$labels) {
        xx <- object$mf[object$id.basis,object$terms[[label]]$vlist]
        x.new <- x[,object$terms[[label]]$vlist]
        nphi <- object$terms[[label]]$nphi
        nrk <-  object$terms[[label]]$nrk
        if (nphi) {
            phi <-  object$terms[[label]]$phi
            for (i in 1:nphi) {
                s <- cbind(s,phi$fun(x.new,nu=i,env=phi$env))
            }
        }
        if (nrk) {
            rk <- object$terms[[label]]$rk
            for (i in 1:nrk) {
                nq <- nq + 1
                r <- r + 10^object$theta[nq]*rk$fun(x.new,xx,nu=i,env=rk$env,out=TRUE)
            }
        }
    }
    as.vector(rho*exp(cbind(s,r)%*%c(object$d,object$c))*object$scal)
}
