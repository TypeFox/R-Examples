## Fit density model
ssden <- function(formula,type=NULL,data=list(),alpha=1.4,
                  weights=NULL,subset,na.action=na.omit,
                  id.basis=NULL,nbasis=NULL,seed=NULL,
                  domain=as.list(NULL),quad=NULL,
                  qdsz.depth=NULL,bias=NULL,
                  prec=1e-7,maxiter=30,skip.iter=FALSE)
{
    ## Obtain model frame and model terms
    mf <- match.call()
    mf$type <- mf$alpha <- NULL
    mf$id.basis <- mf$nbasis <- mf$seed <- NULL
    mf$domain <- mf$quad <- mf$qdsz.depth <- mf$bias <- NULL
    mf$prec <- mf$maxiter <- mf$skip.iter <- NULL
    mf[[1]] <- as.name("model.frame")
    mf <- eval(mf,parent.frame())
    cnt <- model.weights(mf)
    mf$"(weights)" <- NULL
    ## Generate sub-basis
    nobs <- dim(mf)[1]
    if (is.null(id.basis)) {
        if (is.null(nbasis))  nbasis <- max(30,ceiling(10*nobs^(2/9)))
        if (nbasis>=nobs)  nbasis <- nobs
        if (!is.null(seed))  set.seed(seed)
        id.basis <- sample(nobs,nbasis,prob=cnt)
    }
    else {
        if (max(id.basis)>nobs|min(id.basis)<1)
            stop("gss error in ssden: id.basis out of range")
        nbasis <- length(id.basis)
    }
    ## Set domain and/or generate quadrature
    if (is.null(quad)) {
        ## Set domain and type
        fac.list <- NULL
        for (xlab in names(mf)) {
            x <- mf[[xlab]]
            if (is.factor(x)) {
                fac.list <- c(fac.list,xlab)
                domain[[xlab]] <- NULL
            }
            else {
                if (!is.vector(x))
                    stop("gss error in ssden: no default quadrature")
                if (is.null(domain[[xlab]])) {
                    mn <- min(x)
                    mx <- max(x)
                    domain[[xlab]] <- c(mn,mx)+c(-1,1)*(mx-mn)*.05
                }
                else domain[[xlab]] <- c(min(domain[[xlab]]),max(domain[[xlab]]))
                if (is.null(type[[xlab]]))
                    type[[xlab]] <- list("cubic",domain[[xlab]])
                else {
                    if (length(type[[xlab]])==1)
                        type[[xlab]] <- list(type[[xlab]][[1]],domain[[xlab]])
                }
            }
        }
        ## Generate numerical quadrature
        domain <- data.frame(domain)
        mn <- domain[1,]
        mx <- domain[2,]
        dm <- ncol(domain)
        if (dm==1) {
            ## Gauss-Legendre or uniform quadrature
            xlab <- names(domain)
            if (type[[xlab]][[1]]%in%c("per","cubic.per","linear.per")) {
                quad <- list(pt=mn+(1:200)/200*(mx-mn),
                             wt=rep((mx-mn)/200,200))
            }
            else quad <- gauss.quad(200,c(mn,mx))
            quad$pt <- data.frame(quad$pt)
            colnames(quad$pt) <- colnames(domain)
        }
        else {
            ## Smolyak cubature
            if (is.null(qdsz.depth)) qdsz.depth <- switch(min(dm,6)-1,18,14,12,11,10)
            quad <- smolyak.quad(dm,qdsz.depth)
            for (i in 1:ncol(domain)) {
                xlab <- colnames(domain)[i]
                form <- as.formula(paste("~",xlab))
                jk <- ssden(form,data=mf,domain=domain[i],alpha=2,
                            id.basis=id.basis,weights=cnt)
                quad$pt[,i] <- qssden(jk,quad$pt[,i])
                quad$wt <- quad$wt/dssden(jk,quad$pt[,i])
            }
            jk <- NULL
            quad$pt <- data.frame(quad$pt)
            colnames(quad$pt) <- colnames(domain)
        }
        ## Incorporate factors in quadrature
        if (!is.null(fac.list)) {
            for (i in 1:length(fac.list)) {
                wk <-
                  expand.grid(levels(mf[[fac.list[i]]]),1:length(quad$wt))
                quad$wt <- quad$wt[wk[,2]]
                col.names <- c(fac.list[i],colnames(quad$pt))
                quad$pt <- data.frame(wk[,1],quad$pt[wk[,2],])
                colnames(quad$pt) <- col.names
            }
        }
        quad <- list(pt=quad$pt,wt=quad$wt)
    }
    else {
        for (xlab in names(mf)) {
            x <- mf[[xlab]]
            if (is.vector(x)&!is.factor(x)) {
                if (is.null(range <- domain[[xlab]])) {
                    mn <- min(x)
                    mx <- max(x)
                    range <- c(mn,mx)+c(-1,1)*(mx-mn)*.05
                    range[1] <- min(c(range[1],quad$pt[[xlab]]))
                    range[2] <- max(c(range[2],quad$pt[[xlab]]))
                }
                if (is.null(type[[xlab]]))
                    type[[xlab]] <- list("cubic",range)
                else {
                    if (length(type[[xlab]])==1)
                        type[[xlab]] <- list(type[[xlab]][[1]],range)
                    else {
                        mn0 <- min(type[[xlab]][[2]])
                        mx0 <- max(type[[xlab]][[2]])
                        if ((mn0>mn)|(mx0<mx))
                            stop("gss error in ssden: range not covering domain")
                    }
                }
            }
        }
    }
    ## Generate terms
    term <- mkterm(mf,type)
    term$labels <- term$labels[term$labels!="1"]
    ## sampling bias
    qd.pt <- quad$pt
    qd.wt <- quad$wt
    nmesh <- length(qd.wt)
    if (is.null(bias)) {
        nt <- b.wt <- 1
        t.wt <- matrix(1,nmesh,1)
        bias0 <- list(nt=nt,wt=b.wt,qd.wt=t.wt)
    }
    else {
        if ((nt <- length(bias$t))-length(bias$wt))
            stop("gss error in ssden: bias$t and bias$wt mismatch in size")
        b.wt <- abs(bias$wt)/sum(abs(bias$wt))
        t.wt <- NULL
        for (i in 1:nt) t.wt <- cbind(t.wt,bias$fun(bias$t[i],qd.pt))
        bias0 <- list(nt=nt,wt=b.wt,qd.wt=t.wt)
    }
    ## Generate s and r
    s <- qd.s <- r <- qd.r <- NULL
    nq <- 0
    for (label in term$labels) {
        x <- mf[,term[[label]]$vlist]
        x.basis <- mf[id.basis,term[[label]]$vlist]
        qd.x <- qd.pt[,term[[label]]$vlist]
        nphi <- term[[label]]$nphi
        nrk <- term[[label]]$nrk
        if (nphi) {
            phi <- term[[label]]$phi
            for (i in 1:nphi) {
                s <- cbind(s,phi$fun(x,nu=i,env=phi$env))
                qd.s <- cbind(qd.s,phi$fun(qd.x,nu=i,env=phi$env))
            }
        }
        if (nrk) {
            rk <- term[[label]]$rk
            for (i in 1:nrk) {
                nq <- nq+1
                r <- array(c(r,rk$fun(x.basis,x,nu=i,env=rk$env,out=TRUE)),c(nbasis,nobs,nq))
                qd.r <- array(c(qd.r,rk$fun(x.basis,qd.x,nu=i,env=rk$env,out=TRUE)),
                              c(nbasis,nmesh,nq))
            }
        }
    }
    if (!is.null(s)) {
        nnull <- dim(s)[2]
        ## Check s rank
        if (qr(s)$rank<nnull)
            stop("gss error in ssden: unpenalized terms are linearly dependent")
        s <- t(s)
        qd.s <- t(qd.s)
    }
    ## Fit the model
    if (nq==1) {
        r <- r[,,1]
        qd.r <- qd.r[,,1]
        z <- sspdsty(s,r,r[,id.basis],cnt,qd.s,qd.r,qd.wt,prec,maxiter,alpha,bias0)
    }
    else z <- mspdsty(s,r,id.basis,cnt,qd.s,qd.r,qd.wt,prec,maxiter,alpha,bias0,skip.iter)
    ## Brief description of model terms
    desc <- NULL
    for (label in term$labels)
        desc <- rbind(desc,as.numeric(c(term[[label]][c("nphi","nrk")])))
    desc <- rbind(desc,apply(desc,2,sum))
    rownames(desc) <- c(term$labels,"total")
    colnames(desc) <- c("Unpenalized","Penalized")
    ## Return the results
    obj <- c(list(call=match.call(),mf=mf,cnt=cnt,terms=term,desc=desc,
                  alpha=alpha,domain=domain,quad=quad,id.basis=id.basis,
                  qdsz.depth=qdsz.depth,bias=bias0,skip.iter=skip.iter),z)
    class(obj) <- "ssden"
    obj
}

## Fit single smoothing parameter density
sspdsty <- function(s,r,q,cnt,qd.s,qd.r,qd.wt,prec,maxiter,alpha,bias)
{
    nxi <- dim(r)[1]
    nobs <- dim(r)[2]
    nqd <- length(qd.wt)
    if (!is.null(s)) nnull <- dim(s)[1]
    else nnull <- 0
    nxis <- nxi+nnull
    if (is.null(cnt)) cnt <- 0
    ## cv function
    cv <- function(lambda) {
        fit <- .Fortran("dnewton",
                        cd=as.double(cd), as.integer(nxis),
                        as.double(10^(lambda+theta)*q), as.integer(nxi),
                        as.double(rbind(10^theta*r,s)), as.integer(nobs),
                        as.integer(sum(cnt)), as.integer(cnt),
                        as.double(t(rbind(10^theta*qd.r,qd.s))), as.integer(nqd),
                        as.integer(bias$nt), as.double(bias$wt),
                        as.double(t(qd.wt*bias$qd.wt)),
                        as.double(prec), as.integer(maxiter),
                        as.double(.Machine$double.eps), integer(nxis),
                        wk=double(2*((nqd+1)*bias$nt+nobs)+nxis*(2*nxis+4)+max(nxis,3)),
                        info=integer(1),PACKAGE="gss")
        if (fit$info==1) stop("gss error in ssden: Newton iteration diverges")
        if (fit$info==2) warning("gss warning in ssden: Newton iteration fails to converge")
        assign("cd",fit$cd,inherits=TRUE)
        cv <- alpha*fit$wk[2]-fit$wk[1]
        alpha.wk <- max(0,log.la0-lambda-5)*(3-alpha) + alpha
        alpha.wk <- min(alpha.wk,3)
        adj <- ifelse (alpha.wk>alpha,(alpha.wk-alpha)*fit$wk[2],0)
        cv+adj
    }
    ## initialization
    mu.r <- apply(qd.wt*t(qd.r),2,sum)/sum(qd.wt)
    v.r <- apply(qd.wt*t(qd.r^2),2,sum)/sum(qd.wt)
    if (nnull) {
        mu.s <- apply(qd.wt*t(qd.s),2,sum)/sum(qd.wt)
        v.s <- apply(qd.wt*t(qd.s^2),2,sum)/sum(qd.wt)
    }
    if (is.null(s)) theta <- 0
    else theta <- log10(sum(v.s-mu.s^2)/nnull/sum(v.r-mu.r^2)*nxi) / 2
    log.la0 <- log10(sum(v.r-mu.r^2)/sum(diag(q))) + theta
    ## lambda search
    cd <- rep(0,nxi+nnull)
    la <- log.la0
    mn0 <- log.la0-6
    mx0 <- log.la0+6
    repeat {
        mn <- max(la-1,mn0)
        mx <- min(la+1,mx0)
        zz <- nlm0(cv,c(mn,mx))
        if ((min(zz$est-mn,mx-zz$est)>=1e-1)||
            (min(zz$est-mn0,mx0-zz$est)<1e-1)) break
        else la <- zz$est
    }
    ## return
    jk1 <- cv(zz$est)
    int <- sum(qd.wt*exp(t(rbind(10^theta*qd.r,qd.s))%*%cd))
    c <- cd[1:nxi]
    if (nnull) d <- cd[nxi+(1:nnull)]
    else d <- NULL
    list(lambda=zz$est,theta=theta,c=c,d=d,int=int,cv=jk1)
}

## Fit multiple smoothing parameter density
mspdsty <- function(s,r,id.basis,cnt,qd.s,qd.r,qd.wt,prec,maxiter,alpha,bias,skip.iter)
{
    nxi <- dim(r)[1]
    nobs <- dim(r)[2]
    nqd <- length(qd.wt)
    nq <- dim(r)[3]
    if (!is.null(s)) nnull <- dim(s)[1]
    else nnull <- 0
    nxis <- nxi+nnull
    if (is.null(cnt)) cnt <- 0
    ## cv function
    cv <- function(theta) {
        ind.wk <- theta!=theta.old
        if (sum(ind.wk)==nq) {
            r.wk0 <- qd.r.wk0 <- 0
            for (i in 1:nq) {
                r.wk0 <- r.wk0 + 10^theta[i]*r[,,i]
                qd.r.wk0 <- qd.r.wk0 + 10^theta[i]*qd.r[,,i]
            }
            assign("r.wk",r.wk0+0,inherits=TRUE)
            assign("qd.r.wk",qd.r.wk0+0,inherits=TRUE)
            assign("theta.old",theta+0,inherits=TRUE)
        }
        else {
            r.wk0 <- r.wk
            qd.r.wk0 <- qd.r.wk
            for (i in (1:nq)[ind.wk]) {
                theta.wk <- (10^(theta[i]-theta.old[i])-1)*10^theta.old[i]
                r.wk0 <- r.wk0 + theta.wk*r[,,i]
                qd.r.wk0 <- qd.r.wk0 + theta.wk*qd.r[,,i]
            }
        }
        q.wk <- r.wk0[,id.basis]
        fit <- .Fortran("dnewton",
                        cd=as.double(cd), as.integer(nxis),
                        as.double(10^lambda*q.wk), as.integer(nxi),
                        as.double(rbind(r.wk0,s)), as.integer(nobs),
                        as.integer(sum(cnt)), as.integer(cnt),
                        as.double(t(rbind(qd.r.wk0,qd.s))), as.integer(nqd),
                        as.integer(bias$nt), as.double(bias$wt),
                        as.double(t(qd.wt*bias$qd.wt)),
                        as.double(prec), as.integer(maxiter),
                        as.double(.Machine$double.eps), integer(nxis),
                        wk=double(2*((nqd+1)*bias$nt+nobs)+nxis*(2*nxis+4)+max(nxis,3)),
                        info=integer(1),PACKAGE="gss")
        if (fit$info==1) stop("gss error in ssden: Newton iteration diverges")
        if (fit$info==2) warning("gss warning in ssden: Newton iteration fails to converge")
        assign("cd",fit$cd,inherits=TRUE)
        cv <- alpha*fit$wk[2]-fit$wk[1]
        alpha.wk <- max(0,theta-log.th0-5)*(3-alpha) + alpha
        alpha.wk <- min(alpha.wk,3)
        adj <- ifelse (alpha.wk>alpha,(alpha.wk-alpha)*fit$wk[2],0)
        cv+adj
    }
    cv.wk <- function(theta) cv.scale*cv(theta)+cv.shift
    ## initialization
    theta <- -log10(apply(r[,id.basis,],3,function(x)sum(diag(x))))
    r.wk <- qd.r.wk <- 0
    for (i in 1:nq) {
        r.wk <- r.wk + 10^theta[i]*r[,,i]
        qd.r.wk <- qd.r.wk + 10^theta[i]*qd.r[,,i]
    }
    ## theta adjustment
    z <- sspdsty(s,r.wk,r.wk[,id.basis],cnt,qd.s,qd.r.wk,qd.wt,prec,maxiter,alpha,bias)
    theta <- theta + z$theta
    r.wk <- qd.r.wk <- 0
    for (i in 1:nq) {
        theta[i] <- 2*theta[i] + log10(t(z$c)%*%r[,id.basis,i]%*%z$c)
        r.wk <- r.wk + 10^theta[i]*r[,,i]
        qd.r.wk <- qd.r.wk + 10^theta[i]*qd.r[,,i]
    }
    mu <- apply(qd.wt*t(qd.r.wk),2,sum)/sum(qd.wt)
    v <- apply(qd.wt*t(qd.r.wk^2),2,sum)/sum(qd.wt)
    log.la0 <- log10(sum(v-mu^2)/sum(diag(r.wk[,id.basis])))
    log.th0 <- theta-log.la0
    ## lambda search
    z <- sspdsty(s,r.wk,r.wk[,id.basis],cnt,qd.s,qd.r.wk,qd.wt,prec,maxiter,alpha,bias)
    lambda <- z$lambda
    log.th0 <- log.th0 + z$lambda
    theta <- theta + z$theta
    cd <- c(z$c,z$d)
    int <- z$int
    ## early return
    if (skip.iter) {
        z$theta <- theta
        return(z)
    }
    ## theta search
    counter <- 0
    r.wk <- qd.r.wk <- 0
    for (i in 1:nq) {
        r.wk <- r.wk + 10^theta[i]*r[,,i]
        qd.r.wk <- qd.r.wk + 10^theta[i]*qd.r[,,i]
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
        zz <- nlm(cv.wk,theta,stepmax=1,ndigit=7)
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
    qd.r.wk <- 0
    for (i in 1:nq) qd.r.wk <- qd.r.wk + 10^zz$est[i]*qd.r[,,i]
    int <- sum(qd.wt*exp(t(rbind(qd.r.wk,qd.s))%*%cd))
    c <- cd[1:nxi]
    if (nnull) d <- cd[nxi+(1:nnull)]
    else d <- NULL
    list(lambda=lambda,theta=zz$est,c=c,d=d,int=int,cv=jk1)
}
