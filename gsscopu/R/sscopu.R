## Fit copula density model
sscopu <- function(x,symmetry=FALSE,alpha=1.4,order=NULL,exclude=NULL,
                   weights=NULL,id.basis=NULL,nbasis=NULL,seed=NULL,
                   qdsz.depth=NULL,prec=1e-7,maxiter=30,skip.iter=dim(x)[2]!=2)
{
    ## Check inputs
    if ((max(x)>1)|(min(x)<0)) stop("gss error in sscopu: data out of range")
    if (!(is.matrix(x)&dim(x)[2]>=2))
        stop("gss error in sscopu: data must be a matrix of 2 or more columns")
    ## Generate sub-basis
    nobs <- dim(x)[1]
    dm <- dim(x)[2]
    if (is.null(order)) order <- dm
    if (is.null(id.basis)) {
        if (is.null(nbasis))  nbasis <- max(30,ceiling(10*nobs^(2/9)))
        if (nbasis>=nobs)  nbasis <- nobs
        if (!is.null(seed))  set.seed(seed)
        id.basis <- sample(nobs,nbasis,prob=weights)
    }
    else {
        if (max(id.basis)>nobs|min(id.basis)<1)
            stop("gss error in sscopu: id.basis out of range")
        nbasis <- length(id.basis)
    }
    ## exclude
    if (dm==2) exclude <- NULL
    if (!is.null(exclude)) {
        symmetry <- FALSE
        if (is.vector(exclude)) exclude <- matrix(exclude,nrow=1)
        if (dim(exclude)[2]!=2)
            stop("gss error in sscopu: exclude must be a matrix of 2 columns")
        for (i in 1:dim(exclude)[1]) exclude[i,] <- sort(exclude[i,])
        exclude <- unique(exclude)
        if (dim(exclude)[1]==choose(dm,2))
            stop("gss error in sscopu: can not exclude all interactions")
        if (any(exclude[,1]==exclude[,2]))
            stop("gss error in sscopu: only interactions can be excluded")
        if (any(exclude>dm)) stop("gss error in sscopu: exclude out of range")
    }
    ## Generate numerical quadrature
    if (is.null(qdsz.depth)) qdsz.depth <- switch(min(dm,6)-1,24,14,12,11,10)
    quad <- smolyak.quad(dm,qdsz.depth)
    ## Generate terms
    order <- min(dm,max(order,2))
    term <- mkterm.copu(dm,order,symmetry,exclude)
    ## Generate s and r
    s <- qd.s <- r <- qd.r <- NULL
    nq <- 0
    nmesh <- length(quad$wt)
    for (nu in 1:term$nphi) {
        s <- cbind(s,term$phi(x,nu,term$env))
        qd.s <- cbind(qd.s,term$phi(quad$pt,nu,term$env))
    }
    for (nu in 1:term$nrk) {
        nq <- nq+1
        r <- array(c(r,term$rk(x[id.basis,],x,nu,term$env,out=TRUE)),c(nbasis,nobs,nq))
        qd.r <- array(c(qd.r,term$rk(x[id.basis,],quad$pt,nu,term$env,out=TRUE)),
                      c(nbasis,nmesh,nq))
    }
    nnull <- dim(s)[2]
    ## Check s rank
    if (qr(s)$rank<nnull)
        stop("gss error in sscopu: unpenalized terms are linearly dependent")
    s <- t(s)
    qd.s <- t(qd.s)
    ## Fit the model
    nt <- b.wt <- 1
    t.wt <- matrix(1,nmesh,1)
    bias0 <- list(nt=nt,wt=b.wt,qd.wt=t.wt)
    z <- mspdsty(s,r,id.basis,weights,qd.s,qd.r,quad$wt,prec,maxiter,alpha,
                 bias0,skip.iter)
    ## Auxiliary info for copularization
    qd.x <- gauss.quad(200,c(0,1))
    if (dm==2) qd.y <- qd.x
    else {
        qdsz.depth1 <- switch(min(dm-1,6)-1,18,14,12,11,10)
        qd.y <- smolyak.quad(dm-1,qdsz.depth1)
    }
    nn <- 10-1
    x.gd <- (cos((2*(0:nn)+1)/2/(nn+1)*pi)+1)/2
    aa <- cbind(1,2*x.gd-1)
    for (nu in 2:nn) aa <- cbind(aa,(4*x.gd-2)*aa[,nu]-aa[,nu-1])
    if (symmetry) {
        md.wk <- NULL
        for (xx in x.gd) {
            x.wk <- cbind(xx,qd.y$pt)
            s <- NULL
            for (nu in 1:term$nphi) s <- cbind(s,term$phi(x.wk,nu,term$env))
            r <- 0
            for (nu in 1:term$nrk)
                r <- r + 10^z$theta[nu]*term$rk(x.wk,x[id.basis,],nu,term$env,TRUE)
            md.wk <- c(md.wk,sum(exp(s%*%z$d+r%*%z$c)*qd.y$wt)/z$int)
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
        mdsty <- matrix(mdsty,200,dm)
    }
    else {
        mdsty <- NULL
        for (j in 1:dm) {
            md.wk <- NULL
            for (xx in x.gd) {
                x.wk <- matrix(0,length(qd.y$wt),dm)
                x.wk[,j] <- xx
                x.wk[,-j] <- qd.y$pt
                s <- NULL
                for (nu in 1:term$nphi) s <- cbind(s,term$phi(x.wk,nu,term$env))
                r <- 0
                for (nu in 1:term$nrk)
                    r <- r + 10^z$theta[nu]*term$rk(x.wk,x[id.basis,],nu,term$env,TRUE)
                md.wk <- c(md.wk,sum(exp(s%*%z$d+r%*%z$c)*qd.y$wt)/z$int)
            }
            coef <- solve(aa,md.wk)
            tt1 <- 1
            tt2 <- 2*qd.x$pt-1
            mdsty.wk <- coef[1]+coef[2]*tt2
            for (nu in 2:nn) {
                tt.wk <- (4*qd.x$pt-2)*tt2-tt1
                tt1 <- tt2
                tt2 <- tt.wk
                mdsty.wk <- mdsty.wk+coef[nu+1]*tt2
            }
            mdsty.wk <- mdsty.wk/sum(mdsty.wk*qd.x$wt)
            mdsty <- cbind(mdsty,mdsty.wk)
        }
    }
    ## Return the results
    obj <- c(list(call=match.call(),alpha=alpha,id.basis=id.basis,basis=x[id.basis,],
                  env=term$env,nphi=term$nphi,phi=term$phi,nrk=term$nrk,rk=term$rk,
                  mdsty=mdsty,symmetry=symmetry,order=order,qdsz.depth=qdsz.depth),z)
    class(obj) <- "sscopu"
    obj
}
