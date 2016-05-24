## Fit 2-D copula density model, with possibly censored and truncated data
sscopu2 <- function(x,censoring=NULL,truncation=NULL,symmetry=FALSE,alpha=1.4,
                    weights=NULL,id.basis=NULL,nbasis=NULL,seed=NULL,
                    prec=1e-7,maxiter=30)
{
    ## Check inputs
    if ((max(x)>1)|(min(x)<0)) stop("gss error in sscopu2: data out of range")
    if (!(is.matrix(x)&(dim(x)[2]==2)))
        stop("gss error in sscopu2: data must be a matrix of two columns")
    nobs <- dim(x)[1]
    if (!is.null(truncation)) {
        if (!(is.matrix(x)&all(dim(x)==dim(truncation))))
            stop("gss error in sscopu2: truncation and data must match in size")
        if (!all(x<truncation))
            stop("gss error in sscopu2: truncation must be larger than data")
    }
    if (!is.null(censoring)) {
        if (!all(censoring%in%0:3))
            stop("gss error in sscopu2: censoring indicator out of range")
    }
    else censoring <- rep(0,nobs)
    cens <- censoring
    ## Generate sub-basis
    if (is.null(id.basis)) {
        if (is.null(nbasis))  nbasis <- max(30,ceiling(10*nobs^(2/9)))
        if (nbasis>=nobs)  nbasis <- nobs
        if (!is.null(seed))  set.seed(seed)
        id.basis <- sample(nobs,nbasis,prob=weights)
    }
    else {
        if (max(id.basis)>nobs|min(id.basis)<1)
            stop("gss error in sscopu2: id.basis out of range")
        nbasis <- length(id.basis)
    }
    ## Generate numerical quadrature
    hsz <- 20
    qdsz <- 2*hsz
    qd <- gauss.quad(qdsz,c(0,1))
    gap <- diff(qd$pt)
    g.wk <- gap[hsz]/2
    for (i in 1:(hsz-2)) g.wk <- c(g.wk,gap[hsz+i]-g.wk[i])
    g.wk <- 2*g.wk
    g.wk <- c(g.wk,1/2-sum(g.wk))
    gap[hsz:1] <- gap[hsz+(1:hsz)] <- g.wk
    brk <- cumsum(c(0,gap))
    qd.pt <- cbind(rep(qd$pt,qdsz),rep(qd$pt,rep(qdsz,qdsz)))
    nmesh <- qdsz*qdsz
    ## Generate terms
    term <- mkterm.copu(2,2,symmetry,exclude=NULL)
    ## Generate s, r, and q
    idx0 <- (1:length(cens))[cens==0]
    n0 <- length(idx0)
    s0 <- qd.s <- r0 <- q <- qd.r <- NULL
    for (nu in 1:term$nphi) {
        s0 <- cbind(s0,term$phi(x[idx0,],nu,term$env))
        qd.s <- cbind(qd.s,term$phi(qd.pt,nu,term$env))
    }
    nq <- 0
    for (nu in 1:term$nrk) {
        nq <- nq+1
        r0 <- array(c(r0,term$rk(x[id.basis,],x[idx0,],nu,term$env,out=TRUE)),c(nbasis,n0,nq))
        q <- array(c(q,term$rk(x[id.basis,],x[id.basis,],nu,term$env,out=TRUE)),c(nbasis,nbasis,nq))
        qd.r <- array(c(qd.r,term$rk(x[id.basis,],qd.pt,nu,term$env,out=TRUE)),c(nbasis,nmesh,nq))
    }
    idx1 <- (1:length(cens))[cens==1]
    n1 <- length(idx1)
    qd.s1 <- qd.r1 <- wt1 <- NULL
    for (i in idx1) {
        x.wk <- cbind(qd$pt,x[i,2])
        for (nu in 1:term$nphi) qd.s1 <- cbind(qd.s1,term$phi(x.wk,nu,term$env))
        for (nu in 1:term$nrk) qd.r1 <- c(qd.r1,term$rk(x.wk,x[id.basis,],nu,term$env,out=TRUE))
        wt.wk <- qd$wt
        mx <- sum(brk<x[i,1])
        wt.wk[mx:qdsz] <- 0
        wt.wk[mx] <- qd$wt[mx]*(x[i,1]-brk[mx])/gap[mx]
        wt1 <- cbind(wt1,wt.wk)
    }
    if (n1) {
        qd.s1 <- array(qd.s1,c(qdsz,term$nphi,n1))
        qd.r1 <- array(qd.r1,c(qdsz,nbasis,term$nrk,n1))
    }
    idx2 <- (1:length(cens))[cens==2]
    n2 <- length(idx2)
    qd.s2 <- qd.r2 <- wt2 <- NULL
    for (i in idx2) {
        x.wk <- cbind(x[i,1],qd$pt)
        for (nu in 1:term$nphi) qd.s2 <- cbind(qd.s2,term$phi(x.wk,nu,term$env))
        for (nu in 1:term$nrk) qd.r2 <- c(qd.r2,term$rk(x.wk,x[id.basis,],nu,term$env,out=TRUE))
        wt.wk <- qd$wt
        mx <- sum(brk<x[i,2])
        wt.wk[mx:qdsz] <- 0
        wt.wk[mx] <- qd$wt[mx]*(x[i,2]-brk[mx])/gap[mx]
        wt2 <- cbind(wt2,wt.wk)
    }
    if (n2) {
        qd.s2 <- array(qd.s2,c(qdsz,term$nphi,n2))
        qd.r2 <- array(qd.r2,c(qdsz,nbasis,term$nrk,n2))
    }
    idx3 <- (1:length(cens))[cens==3]
    n3 <- length(idx3)
    wt3 <- NULL
    for (i in idx3) {
        for (j in 1:2) {
            wt.wk <- qd$wt
            mx <- sum(brk<x[i,j])
            wt.wk[mx:qdsz] <- 0
            wt.wk[mx] <- qd$wt[mx]*(x[i,j]-brk[mx])/gap[mx]
            wt3 <- cbind(wt3,wt.wk)
        }
    }
    if (n3) wt3 <- array(wt3,c(qdsz,2,n3))
    if (!is.null(weights)) {
        if (n0) cnt0 <- weights[idx0]
        else cnt0 <- NULL
        if (n1) cnt1 <- weights[idx1]
        else cnt1 <- NULL
        if (n2) cnt2 <- weights[idx2]
        else cnt2 <- NULL
        if (n3) cnt3 <- weights[idx3]
        else cnt3 <- NULL
    }
    else {
        cnt0 <- cnt1 <- cnt2 <- cnt3 <- NULL
    }
    ## Group trunccation points
    if (is.null(truncation)) {
        nt <- t.wt <- t.ind <- 1
        wt4 <- array(cbind(qd$wt,qd$wt),c(qdsz,2,nt))
    }
    else {
        wk <- apply(truncation,1,function(x)paste(x,collapse="\r"))
        uwk <- unique(wk)
        nt <- length(uwk)
        wk.dup.ind <- duplicated(wk)
        wk.dup <- as.vector(wk[wk.dup.ind])
        t.ind <- 1:nobs
        t.ind[!wk.dup.ind] <- 1:nt
        if (nobs-nt) {
            t.ind.wk <- range <- 1:(nobs-nt)
            for (i in 1:nt) {
                range.wk <- NULL
                for (j in range) {
                    if (identical(uwk[i],wk.dup[j])) {
                        t.ind.wk[j] <- i
                        range.wk <- c(range.wk,j)
                    }
                }
                if (!is.null(range.wk)) range <- range[!(range%in%range.wk)]
            }
            t.ind[wk.dup.ind] <- t.ind.wk
        }
        t.ind <- t.ind[c(idx0,idx1,idx2,idx3)]
        if (!is.null(weights)) wk <- rep(wk,weights)
        t.wt <- as.vector(table(wk)[uwk])
        t.wt <- t.wt/sum(t.wt)
        wt4 <- NULL
        for (i in (1:nobs)[!wk.dup.ind]) {
            for (j in 1:2) {
                wt.wk <- qd$wt
                mx <- sum(brk<truncation[i,j])
                wt.wk[mx:qdsz] <- 0
                wt.wk[mx] <- qd$wt[mx]*(truncation[i,j]-brk[mx])/gap[mx]
                wt4 <- cbind(wt4,wt.wk)
            }
        }
    }
    trun <- list(nt=nt,t.wt=t.wt,qd.wt=array(wt4,c(qdsz,2,nt)),t.ind=t.ind)
    ## Check s rank
    nnull <- dim(s0)[2]
    if (qr(s0)$rank<nnull)
        stop("gss error in sscopu2: unpenalized terms are linearly dependent")
    s0 <- t(s0)
    qd.s <- t(qd.s)
    ## Fit the model
    z <- mspcopu2(q,s0,r0,cnt0,qd.s,qd.r,qd.s1,qd.r1,wt1,cnt1,qd.s2,qd.r2,wt2,cnt2,
                  wt3,cnt3,trun,prec,maxiter,alpha)
    ## Normalizing constant
    r <- 0
    for (nu in 1:term$nrk) r <- r + 10^z$theta[nu]*t(qd.r[,,nu])
    int <- sum(exp(t(qd.s)%*%z$d+r%*%z$c)*as.vector(outer(qd$wt,qd$wt)))
    ## Auxiliary info for copularization
    qd.x <- gauss.quad(200,c(0,1))
    mdsty <- mdsty2 <- NULL
    nn <- 10-1
    x.gd <- (cos((2*(0:nn)+1)/2/(nn+1)*pi)+1)/2
    aa <- cbind(1,2*x.gd-1)
    for (nu in 2:nn) aa <- cbind(aa,(4*x.gd-2)*aa[,nu]-aa[,nu-1])
    md.wk <- NULL
    for (xx in x.gd) {
        x.wk <- cbind(xx,qd.x$pt)
        s <- NULL
        for (nu in 1:term$nphi) s <- cbind(s,term$phi(x.wk,nu,term$env))
        r <- 0
        for (nu in 1:term$nrk)
            r <- r + 10^z$theta[nu]*term$rk(x.wk,x[id.basis,],nu,term$env,TRUE)
        md.wk <- c(md.wk,sum(exp(s%*%z$d+r%*%z$c)*qd.x$wt)/int)
    }
    coef <- solve(aa,md.wk)
    tt1 <- 1
    tt2 <- 2*qd.x$pt-1
    mdsty <- coef[1]+coef[2]*tt2
    for (nu in 2:nn) {
        tt.wk <- (4*qd.x$pt-2)*tt2-tt1
        tt1 <- tt2
        tt2 <- tt.wk
        mdsty <- mdsty+coef[nu+1]*tt2
    }
    mdsty <- mdsty/sum(mdsty*qd.x$wt)
    if (symmetry) mdsty <- cbind(mdsty,mdsty)
    else {
        md.wk <- NULL
        for (xx in x.gd) {
            x.wk <- cbind(qd.x$pt,xx)
            s <- NULL
            for (nu in 1:term$nphi) s <- cbind(s,term$phi(x.wk,nu,term$env))
            r <- 0
            for (nu in 1:term$nrk)
                r <- r + 10^z$theta[nu]*term$rk(x.wk,x[id.basis,],nu,term$env,TRUE)
            md.wk <- c(md.wk,sum(exp(s%*%z$d+r%*%z$c)*qd.x$wt)/int)
        }
        coef <- solve(aa,md.wk)
        tt1 <- 1
        tt2 <- 2*qd.x$pt-1
        mdsty2 <- coef[1]+coef[2]*tt2
        for (nu in 2:nn) {
            tt.wk <- (4*qd.x$pt-2)*tt2-tt1
            tt1 <- tt2
            tt2 <- tt.wk
            mdsty2 <- mdsty2+coef[nu+1]*tt2
        }
        mdsty2 <- mdsty2/sum(mdsty2*qd.x$wt)
        mdsty <- cbind(mdsty,mdsty2)
    }
    ## Return the results
    obj <- c(list(call=match.call(),alpha=alpha,id.basis=id.basis,basis=x[id.basis,],
                  env=term$env,nphi=term$nphi,phi=term$phi,nrk=term$nrk,rk=term$rk,
                  mdsty=mdsty,symmetry=symmetry,int=int),z)
    class(obj) <- "sscopu"
    obj
}

## Fit multiple smoothing parameter density
mspcopu2 <- function(q,s0,r0,cnt0,qd.s,qd.r,qd.s1,qd.r1,wt1,cnt1,qd.s2,qd.r2,wt2,cnt2,
                     wt3,cnt3,trun,prec,maxiter,alpha)
{
    nxi <- dim(q)[1]
    n0 <- dim(r0)[2]
    nqd <- dim(trun$qd.wt)[1]
    nq <- dim(q)[3]
    nnull <- dim(s0)[1]
    nxis <- nxi+nnull
    if (!is.null(wt1)) n1 <- dim(wt1)[2]
    else n1 <- wt1 <- qd.rs1 <- 0
    if (!is.null(wt2)) n2 <- dim(wt2)[2]
    else n2 <- wt2 <- qd.rs2 <- 0
    if (!is.null(wt3)) n3 <- dim(wt3)[3]
    else n3 <- wt3 <- 0
    if (is.null(cnt0)) cnt0 <- 0
    if (is.null(cnt1)) cnt1 <- 0
    if (is.null(cnt2)) cnt2 <- 0
    if (is.null(cnt3)) cnt3 <- 0
    ## cv function
    cv.s <- function(lambda) {
        fit <- .Fortran("copu2newton",
                        cd=as.double(cd), as.integer(nxis),
                        as.double(10^lambda*q.wk), as.integer(nxi),
                        as.double(rbind(r.wk,s0)), as.integer(n0),
                        as.integer(sum(cnt0)), as.integer(cnt0),
                        as.double(rbind(qd.r.wk,qd.s)), as.integer(nqd),
                        as.double(qd.rs1), as.double(wt1), as.integer(n1),
                        as.integer(sum(cnt1)), as.integer(cnt1),
                        as.double(qd.rs2), as.double(wt2), as.integer(n2),
                        as.integer(sum(cnt2)), as.integer(cnt2),
                        as.double(wt3), as.integer(n3), as.integer(sum(cnt3)),
                        as.integer(cnt3), as.integer(trun$nt), as.double(trun$t.wt),
                        as.double(trun$qd.wt), as.integer(trun$t.ind),
                        as.double(prec), as.integer(maxiter),
                        as.double(.Machine$double.eps), integer(nxis),
                        wk=double(nqd*(2*nqd+n1+n2)+nxis*(2*nxis+trun$nt+5)),
                        info=integer(1),PACKAGE="gsscopu")
        if (fit$info==1) stop("gss error in sscopu2: Newton iteration diverges")
        if (fit$info==2) warning("gss warning in sscopu2: Newton iteration fails to converge")
        assign("cd",fit$cd,inherits=TRUE)
        cv <- alpha*fit$wk[2]+fit$wk[1]
        alpha.wk <- max(0,log.la0-lambda-5)*(3-alpha) + alpha
        alpha.wk <- min(alpha.wk,3)
        adj <- ifelse (alpha.wk>alpha,(alpha.wk-alpha)*fit$wk[2],0)
        cv+adj
    }
    cv <- function(theta) {
        ind.wk <- theta!=theta.old
        if (sum(ind.wk)==nq) {
            q.wk0 <- r.wk0 <- qd.r.wk0 <- 0
            if (n1) qd.r1.wk0 <- 0
            if (n2) qd.r2.wk0 <- 0
            for (i in 1:nq) {
                q.wk0 <- q.wk0 + 10^theta[i]*q[,,i]
                r.wk0 <- r.wk0 + 10^theta[i]*r0[,,i]
                qd.r.wk0 <- qd.r.wk0 + 10^theta[i]*qd.r[,,i]
                if (n1) qd.r1.wk0 <- qd.r1.wk0 + 10^theta[i]*qd.r1[,,i,]
                if (n2) qd.r2.wk0 <- qd.r2.wk0 + 10^theta[i]*qd.r2[,,i,]
            }
            assign("q.wk",q.wk0+0,inherits=TRUE)
            assign("r.wk",r.wk0+0,inherits=TRUE)
            assign("qd.r.wk",qd.r.wk0+0,inherits=TRUE)
            assign("theta.old",theta+0,inherits=TRUE)
            if (n1) assign("qd.r1.wk",qd.r1.wk0+0,inherits=TRUE)
            if (n2) assign("qd.r2.wk",qd.r2.wk0+0,inherits=TRUE)
        }
        else {
            q.wk0 <- q.wk
            r.wk0 <- r.wk
            qd.r.wk0 <- qd.r.wk
            if (n1) qd.r1.wk0 <- qd.r1.wk
            if (n2) qd.r2.wk0 <- qd.r2.wk
            for (i in (1:nq)[ind.wk]) {
                theta.wk <- (10^(theta[i]-theta.old[i])-1)*10^theta.old[i]
                q.wk0 <- q.wk0 + theta.wk*q[,,i]
                r.wk0 <- r.wk0 + theta.wk*r0[,,i]
                qd.r.wk0 <- qd.r.wk0 + theta.wk*qd.r[,,i]
                if (n1) qd.r1.wk0 <- qd.r1.wk0 + theta.wk*qd.r1[,,i,]
                if (n2) qd.r2.wk0 <- qd.r2.wk0 + theta.wk*qd.r2[,,i,]
            }
        }
        if (n1) qd.rs1 <- aperm(array(c(aperm(qd.r1.wk0,c(1,3,2)),
                                        aperm(qd.s1,c(1,3,2))),
                                      c(nqd,n1,nxis)),c(1,3,2))
        if (n2) qd.rs2 <- aperm(array(c(aperm(qd.r2.wk0,c(1,3,2)),
                                        aperm(qd.s2,c(1,3,2))),
                                      c(nqd,n2,nxis)),c(1,3,2))
        fit <- .Fortran("copu2newton",
                        cd=as.double(cd), as.integer(nxis),
                        as.double(10^lambda*q.wk0), as.integer(nxi),
                        as.double(rbind(r.wk0,s0)), as.integer(n0),
                        as.integer(sum(cnt0)), as.integer(cnt0),
                        as.double(rbind(qd.r.wk0,qd.s)), as.integer(nqd),
                        as.double(qd.rs1), as.double(wt1), as.integer(n1),
                        as.integer(sum(cnt1)), as.integer(cnt1),
                        as.double(qd.rs2), as.double(wt2), as.integer(n2),
                        as.integer(sum(cnt2)), as.integer(cnt2),
                        as.double(wt3), as.integer(n3), as.integer(sum(cnt3)),
                        as.integer(cnt3), as.integer(trun$nt), as.double(trun$t.wt),
                        as.double(trun$qd.wt), as.integer(trun$t.ind),
                        as.double(prec), as.integer(maxiter),
                        as.double(.Machine$double.eps), integer(nxis),
                        wk=double(nqd*(2*nqd+n1+n2)+nxis*(2*nxis+trun$nt+5)),
                        info=integer(1),PACKAGE="gsscopu")
        if (fit$info==1) stop("gss error in sscopu2: Newton iteration diverges")
        if (fit$info==2) warning("gss warning in sscopu2: Newton iteration fails to converge")
        assign("cd",fit$cd,inherits=TRUE)
        cv <- alpha*fit$wk[2]+fit$wk[1]
        alpha.wk <- max(0,theta-log.th0-5)*(3-alpha) + alpha
        alpha.wk <- min(alpha.wk,3)
        adj <- ifelse (alpha.wk>alpha,(alpha.wk-alpha)*fit$wk[2],0)
        cv+adj
    }
    cv.wk <- function(theta) cv.scale*cv(theta)+cv.shift
    ## initialization
    theta <- -log10(apply(q,3,function(x)sum(diag(x))))
    q.wk <- r.wk <- qd.r.wk <- 0
    if (n1) qd.r1.wk <- 0
    if (n2) qd.r2.wk <- 0
    for (i in 1:nq) {
        q.wk <- q.wk + 10^theta[i]*q[,,i]
        r.wk <- r.wk + 10^theta[i]*r0[,,i]
        qd.r.wk <- qd.r.wk + 10^theta[i]*qd.r[,,i]
        if (n1) qd.r1.wk <- qd.r1.wk + 10^theta[i]*qd.r1[,,i,]
        if (n2) qd.r2.wk <- qd.r2.wk + 10^theta[i]*qd.r2[,,i,]
    }
    v.r <- v.s <- 0
    for (i in 1:trun$nt) {
        wt.wk <- as.vector(outer(trun$qd.wt[,1,i],trun$qd.wt[,2,i]))
        mu.wk <- apply(t(qd.r.wk)*wt.wk,2,sum)/sum(wt.wk)
        v.r.wk <- apply(t(qd.r.wk^2)*wt.wk,2,sum)/sum(wt.wk)-mu.wk^2
        v.r <- v.r + trun$t.wt[i]*v.r.wk
        mu.wk <- apply(t(qd.s)*wt.wk,2,sum)/sum(wt.wk)
        v.s.wk <- apply(t(qd.s^2)*wt.wk,2,sum)/sum(wt.wk)-mu.wk^2
        v.s <- v.s + trun$t.wt[i]*v.s.wk
    }
    log.la0 <- log10(sum(v.r)/sum(diag(q.wk)))
    ## initial lambda search
    theta.wk <- log10(sum(v.s)/nnull/sum(v.r)*nxi) / 2
    q.wk <- 10^theta.wk*q.wk
    r.wk <- 10^theta.wk*r.wk
    qd.r.wk <- 10^theta.wk*qd.r.wk
    if (n1) qd.r1.wk <- 10^theta.wk*qd.r1.wk
    if (n2) qd.r2.wk <- 10^theta.wk*qd.r2.wk
    theta <- theta + theta.wk
    log.la0 <- log.la0 + theta.wk
    if (n1) qd.rs1 <- aperm(array(c(aperm(qd.r1.wk,c(1,3,2)),aperm(qd.s1,c(1,3,2))),
                                  c(nqd,n1,nxis)),c(1,3,2))
    if (n2) qd.rs2 <- aperm(array(c(aperm(qd.r2.wk,c(1,3,2)),aperm(qd.s2,c(1,3,2))),
                                  c(nqd,n2,nxis)),c(1,3,2))
    cd <- rep(0,nxi+nnull)
    la <- log.la0
    repeat {
        mn <- la-1
        mx <- la+1
        if (mx>log.la0+6) break
        zz <- nlm0(cv.s,c(mn,mx))
        if (min(zz$est-mn,mx-zz$est)>=1e-3) break
        else la <- zz$est
    }
    ## theta adjustment
    q.wk <- r.wk <- qd.r.wk <- 0
    if (n1) qd.r1.wk <- 0
    if (n2) qd.r2.wk <- 0
    for (i in 1:nq) {
        theta[i] <- 2*theta[i] + log10(t(cd[1:nxi])%*%q[,,i]%*%cd[1:nxi])
        q.wk <- q.wk + 10^theta[i]*q[,,i]
        r.wk <- r.wk + 10^theta[i]*r0[,,i]
        qd.r.wk <- qd.r.wk + 10^theta[i]*qd.r[,,i]
        if (n1) qd.r1.wk <- qd.r1.wk + 10^theta[i]*qd.r1[,,i,]
        if (n2) qd.r2.wk <- qd.r2.wk + 10^theta[i]*qd.r2[,,i,]
    }
    v.r <- v.s <- 0
    for (i in 1:trun$nt) {
        wt.wk <- as.vector(outer(trun$qd.wt[,1,i],trun$qd.wt[,2,i]))
        if (sum(wt.wk)==0) next
        mu.wk <- apply(t(qd.r.wk)*wt.wk,2,sum)/sum(wt.wk)
        v.r.wk <- apply(t(qd.r.wk^2)*wt.wk,2,sum)/sum(wt.wk)-mu.wk^2
        v.r <- v.r + trun$t.wt[i]*v.r.wk
        mu.wk <- apply(t(qd.s)*wt.wk,2,sum)/sum(wt.wk)
        v.s.wk <- apply(t(qd.s^2)*wt.wk,2,sum)/sum(wt.wk)-mu.wk^2
        v.s <- v.s + trun$t.wt[i]*v.s.wk
    }
    log.la0 <- log10(sum(v.r)/sum(diag(q.wk)))
    ## lambda search
    theta.wk <- log10(sum(v.s)/nnull/sum(v.r)*nxi) / 2
    q.wk <- 10^theta.wk*q.wk
    r.wk <- 10^theta.wk*r.wk
    qd.r.wk <- 10^theta.wk*qd.r.wk
    if (n1) qd.r1.wk <- 10^theta.wk*qd.r1.wk
    if (n2) qd.r2.wk <- 10^theta.wk*qd.r2.wk
    theta <- theta + theta.wk
    log.la0 <- log.la0 + theta.wk
    if (n1) qd.rs1 <- aperm(array(c(aperm(qd.r1.wk,c(1,3,2)),aperm(qd.s1,c(1,3,2))),
                                  c(nqd,n1,nxis)),c(1,3,2))
    if (n2) qd.rs2 <- aperm(array(c(aperm(qd.r2.wk,c(1,3,2)),aperm(qd.s2,c(1,3,2))),
                                  c(nqd,n2,nxis)),c(1,3,2))
    cd <- rep(0,nxi+nnull)
    la <- log.la0
    repeat {
        mn <- la-1
        mx <- la+1
        if (mx>log.la0+6) break
        zz <- nlm0(cv.s,c(mn,mx))
        if (min(zz$est-mn,mx-zz$est)>=1e-3) break
        else la <- zz$est
    }
    log.th0 <- theta-log.la0
    lambda <- zz$est
    log.th0 <- log.th0 + lambda
    ## theta search
    counter <- 0
    q.wk <- r.wk <- qd.r.wk <- 0
    if (n1) qd.r1.wk <- 0
    if (n2) qd.r2.wk <- 0
    for (i in 1:nq) {
        q.wk <- q.wk + 10^theta[i]*q[,,i]
        r.wk <- r.wk + 10^theta[i]*r0[,,i]
        qd.r.wk <- qd.r.wk + 10^theta[i]*qd.r[,,i]
        if (n1) qd.r1.wk <- qd.r1.wk + 10^theta[i]*qd.r1[,,i,]
        if (n2) qd.r2.wk <- qd.r2.wk + 10^theta[i]*qd.r2[,,i,]
    }
    theta.old <- theta
    ## scale and shift cv
    tmp <- abs(cv(theta))
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
    repeat {
        zz <- nlm(cv.wk,theta,stepmax=.5,ndigit=7)
        if (zz$code<=3)  break
        theta <- zz$est        
        counter <- counter + 1
        if (counter>=5) {
            warning("gss warning in ssden: CV iteration fails to converge")
            break
        }
    }
    ## return
    jk1 <- cv(zz$est)
    c <- cd[1:nxi]
    d <- cd[nxi+(1:nnull)]
    list(lambda=lambda,theta=zz$est,c=c,d=d,cv=jk1)
}
