dsscopu <-  ## Evaluate copula density
function(object,x,copu=TRUE) {
    ## Check inputs
    if (class(object)!="sscopu")
        stop("gss error in dsscopu: not a sscopu object")
    if ((max(x)>1)|(min(x)<0)) stop("gss error in dsscopu: points out of range")
    if (is.vector(x)) x <- matrix(x,1,dim(object$basis)[2])
    ## obtain quantiles of marginal distribution
    if (copu) {
        md <- object$mdsty
        qd <- gauss.quad(200,c(0,1))
        gap <- diff(qd$pt)
        g.wk <- gap[100]/2
        for (i in 1:98) g.wk <- c(g.wk,gap[100+i]-g.wk[i])
        g.wk <- 2*g.wk
        g.wk <- c(g.wk,1/2-sum(g.wk))
        gap[100:1] <- gap[100+(1:100)] <- g.wk
        brk <- cumsum(c(0,gap))
        qq <- NULL
        for (j in 1:dim(md)[2]) {
            p <- x[,j]
            order.p <- rank(p)
            q <- p <- sort(p)
            p.dup <- duplicated(p)
            p.wk <- cumsum(md[,j]*qd$wt)
            kk <- (1:length(p))[p>0&p<1]
            for (i in kk) {
                if (p.dup[i]) {
                    q[i] <- q.dup
                    next
                }
                ind <- (1:200)[p.wk<p[i]]
                if (!length(ind)) wk <- p[i]/(md[1,j]*qd$wt[1])*gap[1]
                else {
                    id.mx <- max(ind)
                    wk <- brk[id.mx+1]+(p[i]-p.wk[id.mx])/
                      (md[id.mx+1,j]*qd$wt[id.mx+1])*gap[id.mx+1]
                }
                q[i] <- q.dup <- wk
            }
            qq <- cbind(qq,q[order.p])
        }
    }
    else qq <- x
    ## evaluate density
    s <- NULL
    for (nu in 1:object$nphi) s <- cbind(s,object$phi(qq,nu,object$env))
    r <- 0
    for (nu in 1:object$nrk)
        r <- r + 10^object$theta[nu]*object$rk(qq,object$basis,nu,object$env,TRUE)
    z <- as.vector(exp(s%*%object$d+r%*%object$c))/object$int
    if (copu) {
        for (i in 1:dim(qq)[1]) {
            for (j in 1:dim(md)[2]) {
                ind <- order(abs(qd$pt-qq[i,j]))[1:4]
                dd <- qd$pt[ind]-qq[i,j]
                aa <- cbind(1,dd,dd^2,dd^3)
                wk <- solve(aa,md[ind,j])[1]
                z[i] <- z[i]/wk
            }
        }
    }
    z
}

cdsscopu <- ## Evaluate 1-D conditional density
function (object,x,cond,pos=1,int=NULL) {
    if (class(object)!="sscopu")
        stop("gss error in cdsscopu: not a sscopu object")
    dm <- dim(object$basis)[2]
    if (length(cond)!=dm-1)
        stop("gss error in cdsscopu: condition is of wrong dimension")
    if ((min(x,cond)<0)|(max(x,cond)>1))
        stop("gss error in cdsscopu: points out of range")
    if ((pos<1)|(pos>dm)) stop("gss error in cdsscopu: position out of range")
    if (is.null(int)) {
        quad <- gauss.quad(200,c(0,1))
        xx <- matrix(0,200,dm)
        xx[,-pos] <- t(matrix(cond,dm-1,200))
        xx[,pos] <- quad$pt
        int <- sum(dsscopu(object,xx)*quad$wt)
    }
    ## Return value
    xx <- matrix(0,length(x),dm)
    xx[,-pos] <- t(matrix(cond,dm-1,length(x)))
    xx[,pos] <- x
    list(pdf=dsscopu(object,xx)/int,int=int)
}

cpsscopu <- ## Compute cdf for 1-D conditional density
function(object,q,cond,pos=1) {
    if (class(object)!="sscopu")
        stop("gss error in cpsscopu: not a sscopu object")
    dm <- dim(object$basis)[2]
    if (length(cond)!=dm-1)
        stop("gss error in cpsscopu: condition is of wrong dimension")
    if ((min(q,cond)<0)|(max(q,cond)>1))
        stop("gss error in cpsscopu: points out of range")
    if ((pos<1)|(pos>dm)) stop("gss error in cpsscopu: position out of range")
    order.q <- rank(q)
    p <- q <- sort(q)
    q.dup <- duplicated(q)
    p[q<=0] <- 0
    p[q>=1] <- 1
    qd.hize <- 200
    qd <- gauss.quad(2*qd.hize,c(0,1))
    d.qd <- cdsscopu(object,qd$pt,cond,pos)$pdf
    gap <- diff(qd$pt)
    g.wk <- gap[qd.hize]/2
    for (i in 1:(qd.hize-2)) g.wk <- c(g.wk,gap[qd.hize+i]-g.wk[i])
    g.wk <- 2*g.wk
    g.wk <- c(g.wk,1/2-sum(g.wk))
    gap[qd.hize:1] <- gap[qd.hize+(1:qd.hize)] <- g.wk
    brk <- cumsum(c(0,gap))
    kk <- (1:length(q))[q>0&q<1]
    for (i in kk) {
        if (q.dup[i]) {
            p[i] <- p.dup
            next
        }
        ind <- (1:(2*qd.hize))[qd$pt<q[i]]
        if (!length(ind)) {
            wk <- d.qd[1]*qd$wt[1]*q[i]/gap[1]
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

cqsscopu <- ## Compute quantiles for 1-D conditional density
function(object,p,cond,pos=1) {
    if (class(object)!="sscopu")
        stop("gss error in cqsscopu: not a sscopu object")
    dm <- dim(object$basis)[2]
    if (length(cond)!=dm-1)
        stop("gss error in cqsscopu: condition is of wrong dimension")
    if ((min(p,cond)<0)|(max(p,cond)>1))
        stop("gss error in cqsscopu: points out of range")
    if ((pos<1)|(pos>dm)) stop("gss error in cqsscopu: position out of range")
    order.p <- rank(p)
    q <- p <- sort(p)
    p.dup <- duplicated(p)
    q[p<=0] <- 0
    q[p>=1] <- 1
    qd.hize <- 200
    qd <- gauss.quad(2*qd.hize,c(0,1))
    d.qd <- cdsscopu(object,qd$pt,cond,pos)$pdf
    gap <- diff(qd$pt)
    g.wk <- gap[qd.hize]/2
    for (i in 1:(qd.hize-2)) g.wk <- c(g.wk,gap[qd.hize+i]-g.wk[i])
    g.wk <- 2*g.wk
    g.wk <- c(g.wk,1/2-sum(g.wk))
    gap[qd.hize:1] <- gap[qd.hize+(1:qd.hize)] <- g.wk
    brk <- cumsum(c(0,gap))
    p.wk <- cumsum(d.qd*qd$wt)
    kk <- (1:length(p))[p>0&p<1]
    for (i in kk) {
        if (p.dup[i]) {
            q[i] <- q.dup
            next
        }
        ind <- (1:(2*qd.hize))[p.wk<p[i]]
        if (!length(ind)) {
            wk <- p[i]/(d.qd[1]*qd$wt[1])*gap[1]
        }
        else {
            id.mx <- max(ind)
            wk <- brk[id.mx+1]+(p[i]-p.wk[id.mx])/(d.qd[id.mx+1]*qd$wt[id.mx+1])*gap[id.mx+1]
        }
        q[i] <- q.dup <- wk
    }
    q[order.p]
}
