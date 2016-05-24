cdssden <- ## Evaluate conditional density
function (object,x,cond,int=NULL) {
    if (!("ssden"%in%class(object))) stop("gss error in cdssden: not a ssden object")
    if (nrow(cond)!=1) stop("gss error in cdssden: condition has to be a single point")
    xnames <- NULL
    for (i in colnames(object$mf))
        if (all(i!=colnames(cond))) xnames <- c(xnames,i)
    if (any(length(xnames)==c(0,ncol(object$mf))))
        stop("gss error in cdssden: not a conditional density")
    if (length(xnames)==1&is.vector(x)) {
        x <- data.frame(x)
        colnames(x) <- xnames
    }
    if (!all(sort(xnames)==sort(colnames(x))))
        stop("gss error in cdssden: mismatched variable names")
    ## Calculate normalizing constant
    while (is.null(int)) {
        fac.list <- NULL
        num.list <- NULL
        for (xlab in xnames) {
            if (is.factor(x.wk <- x[[xlab]])) fac.list <- c(fac.list,xlab)
            else {
                if (!is.vector(x.wk)|is.null(object$domain[[xlab]])) {
                    warning("gss warning in cdssden: int set to 1")
                    int <- 1
                    next
                }
                else num.list <- c(num.list,xlab)
            }
        }
        ## Generate quadrature for numerical variables
        if (!is.null(num.list)) {
            if (length(num.list)==1) {
                ## Gauss-Legendre quadrature
                mn <- min(object$domain[,num.list])
                mx <- max(object$domain[,num.list])
                quad <- gauss.quad(200,c(mn,mx))
                quad$pt <- data.frame(quad$pt)
                colnames(quad$pt) <- num.list
            }
            else {
                ## Smolyak cubature
                domain.wk <- object$domain[,num.list]
                code <- c(15,14,13)
                quad <- smolyak.quad(ncol(domain.wk),code[ncol(domain.wk)-1])
                for (i in 1:ncol(domain.wk)) {
                    xlab <- colnames(domain.wk)[i]
                    form <- as.formula(paste("~",xlab))
                    jk <- ssden(form,data=object$mf,domain=domain.wk[i],alpha=2,
                                id.basis=object$id.basis)
                    quad$pt[,i] <- qssden(jk,quad$pt[,i])
                    quad$wt <- quad$wt/dssden(jk,quad$pt[,i])
                }
                jk <- NULL
                quad$pt <- data.frame(quad$pt)
                colnames(quad$pt) <- colnames(domain.wk)
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
        xmesh <- quad$pt[,xnames,drop=FALSE]
        xx <- cond[rep(1,nrow(xmesh)),,drop=FALSE]
        int <- sum(dssden(object,cbind(xmesh,xx))*quad$wt)
    }
    ## Return value
    xx <- cond[rep(1,nrow(x)),,drop=FALSE]
    list(pdf=dssden(object,cbind(x,xx))/int,int=int)
}

cpssden <- ## Compute cdf for univariate conditional density
function(object,q,cond) {
    if (!("ssden"%in%class(object))) stop("gss error in cpssden: not a ssden object")
    xnames <- NULL
    for (i in colnames(object$mf))
        if (all(i!=colnames(cond))) xnames <- c(xnames,i)
    if ((length(xnames)!=1)|!is.vector(object$mf[,xnames]))
        stop("gss error in cpssden: not a 1-D conditional density")
    mn <- min(object$domain[,xnames])
    mx <- max(object$domain[,xnames])
    order.q <- rank(q)
    p <- q <- sort(q)
    q.dup <- duplicated(q)
    p[q<=mn] <- 0
    p[q>=mx] <- 1
    qd.hize <- 200
    qd <- gauss.quad(2*qd.hize,c(mn,mx))
    d.qd <- cdssden(object,qd$pt,cond)$pdf
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

cqssden <- ## Compute quantiles for univariate conditional density
function(object,p,cond) {
    if (!("ssden"%in%class(object))) stop("gss error in cqssden: not a ssden object")
    xnames <- NULL
    for (i in colnames(object$mf))
        if (all(i!=colnames(cond))) xnames <- c(xnames,i)
    if ((length(xnames)!=1)|!is.vector(object$mf[,xnames]))
        stop("gss error in cqssden: not a 1-D conditional density")
    mn <- min(object$domain[,xnames])
    mx <- max(object$domain[,xnames])
    order.p <- rank(p)
    q <- p <- sort(p)
    p.dup <- duplicated(p)
    q[p<=0] <- mn
    q[p>=1] <- mx
    qd.hize <- 200
    qd <- gauss.quad(2*qd.hize,c(mn,mx))
    d.qd <- cdssden(object,qd$pt,cond)$pdf
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
