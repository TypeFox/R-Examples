cdsscden <- ## Evaluate conditional density estimate
function (object,y,x,cond,int=NULL) {
    ## check inputs
    if (!("sscden"%in%class(object))) stop("gss error in cdsscden: not a sscden object")
    if (!all(sort(object$xnames)==sort(colnames(x))))
        stop("gss error in cdsscden: mismatched x variable names")
    if (nrow(cond)!=1) stop("gss error in cdsscden: condition has to be a single point")
    ynames <- NULL
    for (i in object$ynames)
        if (all(i!=colnames(cond))) ynames <- c(ynames,i)
    if (any(length(ynames)==c(0,length(object$ynames))))
        stop("gss error in cdsscden: not a conditional density")
    if (length(ynames)==1&is.vector(y)) {
        y <- data.frame(y)
        colnames(y) <- ynames
    }
    if (!all(sort(ynames)==sort(colnames(y))))
        stop("gss error in cdsscden: mismatched y variable names")
    ## Calculate normalizing constant
    if ("sscden1"%in%class(object)) ydomain <- object$rho$env$ydomain
    else ydomain <- object$ydomain
    while (is.null(int)) {
        fac.list <- NULL
        num.list <- NULL
        for (ylab in ynames) {
            if (is.factor(y.wk <- y[[ylab]])) fac.list <- c(fac.list,ylab)
            else {
                if (!is.vector(y.wk)|is.null(ydomain[[ylab]])) {
                    warning("gss warning in cdsscden: int set to 1")
                    int <- 1
                    next
                }
                else num.list <- c(num.list,ylab)
            }
        }
        ## Generate quadrature for numerical variables
        if (!is.null(num.list)) {
            if (length(num.list)==1) {
                ## Gauss-Legendre quadrature
                mn <- min(ydomain[,num.list])
                mx <- max(ydomain[,num.list])
                quad <- gauss.quad(200,c(mn,mx))
                quad$pt <- data.frame(quad$pt)
                colnames(quad$pt) <- num.list
            }
            else {
                ## Smolyak cubature
                domain.wk <- ydomain[,num.list]
                code <- c(15,14,13)
                quad <- smolyak.quad(ncol(domain.wk),code[ncol(domain.wk)-1])
                for (i in 1:ncol(domain.wk)) {
                    ylab <- colnames(domain.wk)[i]
                    wk <- object$mf[[ylab]]
                    jk <- ssden(~wk,domain=data.frame(wk=domain.wk[,i]),alpha=2,
                                id.basis=object$id.basis)
                    quad$pt[,i] <- qssden(jk,quad$pt[,i])
                    quad$wt <- quad$wt/dssden(jk,quad$pt[,i])
                }
                jk <- wk <- NULL
                quad$pt <- data.frame(quad$pt)
                colnames(quad.pt) <- colnames(domain.wk)
            }
        }
        else quad <- list(pt=data.frame(dum=1),wt=1)
        ## Incorporate factors in quadrature
        if (!is.null(fac.list)) {
            for (i in 1:length(fac.list)) {
                wk <- expand.grid(levels(object$mf[[fac.list[i]]]),1:length(quad$wt))
                quad$wt <- quad$wt[wk[,2]]
                col.names <- c(fac.list[i],colnames(quad$pt))
                quad$pt <- data.frame(wk[,1],quad$pt[wk[,2],])
                colnames(quad$pt) <- col.names
            }
        }
        ymesh <- quad$pt[,ynames,drop=FALSE]
        yy <- cond[rep(1,nrow(ymesh)),,drop=FALSE]
        int <- apply(dsscden(object,cbind(ymesh,yy),x)*quad$wt,2,sum)
    }
    ## Return value
    yy <- cond[rep(1,nrow(y)),,drop=FALSE]
    list(pdf=t(t(dsscden(object,cbind(y,yy),x))/int),int=int)
}

cpsscden <- ## Compute cdf for univariate conditional density
function(object,q,x,cond) {
    ## check inputs
    if (!("sscden"%in%class(object))) stop("gss error in cpsscden: not a sscden object")
    if (!all(sort(object$xnames)==sort(colnames(x))))
        stop("gss error in cpsscden: mismatched x variable names")
    if (nrow(cond)!=1) stop("gss error in cpsscden: condition has to be a single point")
    ynames <- NULL
    for (i in object$ynames) if (all(i!=colnames(cond))) ynames <- c(ynames,i)
    if (length(ynames)!=1) stop("gss error in cpsscden: y is not 1-D")
    if (is.factor(object$mf[,ynames])) stop("gss error in cpsscden: y is not continuous")
    if ("sscden1"%in%class(object)) ydomain <- object$rho$env$ydomain
    else ydomain <- object$ydomain
    mn <- min(ydomain[[ynames]])
    mx <- max(ydomain[[ynames]])
    order.q <- rank(q)
    p <- q <- sort(q)
    q.dup <- duplicated(q)
    p[q<=mn] <- 0
    p[q>=mx] <- 1
    qd.hize <- 200
    qd <- gauss.quad(2*qd.hize,c(mn,mx))
    y.wk <- data.frame(qd$pt)
    colnames(y.wk) <- ynames
    y.wk <- cbind(y.wk,cond)
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

cqsscden <- ## Compute quantiles for univariate density estimate
function(object,p,x,cond) {
    ## check inputs
    if (!("sscden"%in%class(object))) stop("gss error in cqsscden: not a sscden object")
    if (!all(sort(object$xnames)==sort(colnames(x))))
        stop("gss error in cqsscden: mismatched x variable names")
    if (nrow(cond)!=1) stop("gss error in cqsscden: condition has to be a single point")
    ynames <- NULL
    for (i in object$ynames) if (all(i!=colnames(cond))) ynames <- c(ynames,i)
    if (length(ynames)!=1) stop("gss error in cqsscden: y is not 1-D")
    if (is.factor(object$mf[,ynames])) stop("gss error in cqsscden: y is not continuous")
    if ("sscden1"%in%class(object)) ydomain <- object$rho$env$ydomain
    else ydomain <- object$ydomain
    mn <- min(ydomain[[ynames]])
    mx <- max(ydomain[[ynames]])
    order.p <- rank(p)
    q <- p <- sort(p)
    p.dup <- duplicated(p)
    q[p<=0] <- mn
    q[p>=1] <- mx
    qd.hize <- 200
    qd <- gauss.quad(2*qd.hize,c(mn,mx))
    y.wk <- data.frame(qd$pt)
    colnames(y.wk) <- ynames
    y.wk <- cbind(y.wk,cond)
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
