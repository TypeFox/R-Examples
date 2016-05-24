## Fit log-linear regression model
sscden <- function(formula,response,type=NULL,data=list(),weights,
                   subset,na.action=na.omit,alpha=1.4,
                   id.basis=NULL,nbasis=NULL,seed=NULL,
                   ydomain=as.list(NULL),yquad=NULL,
                   prec=1e-7,maxiter=30,skip.iter=FALSE)
{
    ## Obtain model frame and model terms
    mf <- match.call()
    mf$response <- mf$type <- mf$alpha <- NULL
    mf$id.basis <- mf$nbasis <- mf$seed <- NULL
    mf$ydomain <- mf$yquad <- NULL
    mf$prec <- mf$maxiter <- mf$skip.iter <- NULL
    term.wk <- terms.formula(formula)
    ynames <- as.character(attr(terms(response),"variables"))[-1]
    mf[[1]] <- as.name("model.frame")
    mf <- eval(mf,parent.frame())
    nobs <- nrow(mf)
    cnt <- model.weights(mf)
    if (is.null(cnt)) data$cnt <- rep(1,nobs)
    else {
        data$cnt <- cnt
        mf$"(weights)" <- NULL
    }
    ## Generate sub-basis
    nobs <- nrow(mf)
    if (is.null(id.basis)) {
        if (is.null(nbasis))  nbasis <- max(30,ceiling(10*nobs^(2/9)))
        if (nbasis>=nobs)  nbasis <- nobs
        if (!is.null(seed))  set.seed(seed)
        id.basis <- sample(nobs,nbasis,prob=cnt)
    }
    else {
        if (max(id.basis)>nobs|min(id.basis)<1)
            stop("gss error in sscden: id.basis out of range")
        nbasis <- length(id.basis)
    }
    ## Check inputs
    mt <- attr(mf,"terms")
    vars <- as.character(attr(mt,"variables"))[-1]
    if(!all(ynames%in%vars)) stop("gss error in sscden: response missing in model")
    xnames <- vars[!(vars%in%ynames)]
    if (is.null(xnames)) stop("gss error in sscden: missing covariate")
    ## Set ydomain and type
    mtrx.y <- FALSE
    for (ylab in ynames) {
        y <- mf[[ylab]]
        if (!is.factor(y)) {
            if (is.vector(y)) {
                if (is.null(ydomain[[ylab]])) {
                    mn <- min(y)
                    mx <- max(y)
                    ydomain[[ylab]] <- c(mn,mx)+c(-1,1)*(mx-mn)*.05
                }
                else ydomain[[ylab]] <- c(min(ydomain[[ylab]]),max(ydomain[[ylab]]))
                if (is.null(type[[ylab]]))
                    type[[ylab]] <- list("cubic",ydomain[[ylab]])
                else {
                    if (length(type[[ylab]])==1)
                        type[[ylab]] <- list(type[[ylab]][[1]],ydomain[[ylab]])
                }
            }
            else mtrx.y <- TRUE
        }
    }
    ydomain <- data.frame(ydomain)
    ## Generate terms    
    term <- mkterm(mf,type)
    term.labels <- labels(mt)
    facs <- attr(mt,"factors")
    ind.wk <- NULL
    for (lab in term.labels)
        ind.wk <- c(ind.wk,any(facs[ynames,lab]))
    term$labels <- term.labels[ind.wk]
    ## Generate quadrature
    if (is.null(yquad)) {
        if (mtrx.y) stop("gss error in sscden: no default quadrature")
        yquad <- ssden(response,id.basis=id.basis,data=data,weights=cnt,
                       alpha=2,domain=ydomain)$quad
    }
    qd.pt <- yquad$pt
    qd.wt <- yquad$wt
    nmesh <- length(qd.wt)
    ## obtain unique covariate observations
    x <- xx <- mf[,xnames,drop=FALSE]
    xx <- apply(xx,1,function(x)paste(x,collapse="\r"))
    x.dup.ind <- duplicated(xx)    
    if (!is.null(cnt)) xx <- rep(xx,cnt)
    xx.wt <- as.vector(table(xx)[unique(xx)])
    xx.wt <- xx.wt/sum(xx.wt)
    nx <- length(xx.wt)
    ## Generate s, r, qd.s, and qd.r
    s <- r <- qd.s <- NULL
    qd.r <- as.list(NULL)
    nu <- nq <- 0
    for (label in term$labels) {
        vlist <- term[[label]]$vlist
        x.list <- xnames[xnames%in%vlist]
        y.list <- ynames[ynames%in%vlist]
        xy <- mf[,vlist]
        xy.basis <- mf[id.basis,vlist]
        qd.xy <- data.frame(matrix(0,nmesh,length(vlist)))
        names(qd.xy) <- vlist
        qd.xy[,y.list] <- qd.pt[,y.list]
        if (length(x.list)) xx <- x[!x.dup.ind,x.list,drop=FALSE]
        else xx <- NULL
        nphi <- term[[label]]$nphi
        nrk <- term[[label]]$nrk
        if (nphi) {
            phi <- term[[label]]$phi
            for (i in 1:nphi) {
                nu <- nu+1
                s.wk <- phi$fun(xy,nu=i,env=phi$env)
                s <- cbind(s,s.wk)
                if (is.null(xx)) {
                    qd.s.wk <- phi$fun(qd.xy[,,drop=TRUE],nu=i,env=phi$env)
                    qd.wk <- matrix(qd.s.wk,nmesh,nx)
                }
                else {
                    qd.wk <- NULL
                    for (j in 1:nx) {
                        qd.xy[,x.list] <- xx[rep(j,nmesh),]
                        qd.wk <- cbind(qd.wk,phi$fun(qd.xy,i,phi$env))
                    }
                }
                qd.s <- array(c(qd.s,qd.wk),c(nmesh,nx,nu))
            }
        }
        if (nrk) {
            rk <- term[[label]]$rk
            for (i in 1:nrk) {
                nq <- nq+1
                r.wk <- rk$fun(xy,xy.basis,nu=i,env=rk$env,out=TRUE)
                r <- array(c(r,r.wk),c(nobs,nbasis,nq))
                if (is.null(xx)) {
                    qd.r.wk <- rk$fun(qd.xy[,,drop=TRUE],xy.basis,nu=i,env=rk$env,out=TRUE)
                    qd.r[[nq]] <- qd.r.wk
                }
                else {
                    qd.wk <- NULL
                    for (j in 1:nx) {
                        qd.xy[,x.list] <- xx[rep(j,nmesh),]
                        qd.wk <- array(c(qd.wk,rk$fun(qd.xy,xy.basis,i,rk$env,TRUE)),
                                       c(nmesh,nbasis,j))
                    }
                    qd.r[[nq]] <- qd.wk                    
                }
            }
        }
    }
    ## Check s rank
    if (!is.null(s)) {
        nnull <- dim(s)[2]
        if (qr(s)$rank<nnull)
            stop("gss error in sscden: unpenalized MLE is not unique")
    }
    ## Fit the model
    z <- mspcdsty(s,r,id.basis,cnt,qd.s,qd.r,xx.wt,qd.wt,prec,maxiter,alpha,skip.iter)
    ## Brief description of model terms
    desc <- NULL
    for (label in term$labels)
        desc <- rbind(desc,as.numeric(c(term[[label]][c("nphi","nrk")])))
    desc <- rbind(desc,apply(desc,2,sum))
    rownames(desc) <- c(term$labels,"total")
    colnames(desc) <- c("Unpenalized","Penalized")
    ## Return the results
    obj <- c(list(call=match.call(),mf=mf,cnt=cnt,terms=term,desc=desc,
                  ydomain=ydomain,yquad=yquad,xx.wt=xx.wt,x.dup.ind=x.dup.ind,
                  alpha=alpha,ynames=ynames,xnames=xnames,id.basis=id.basis,
                  skip.iter=skip.iter),z)
    class(obj) <- c("sscden")
    obj
}

## Fit (multiple smoothing parameter) log-linear regression model
mspcdsty <- function(s,r,id.basis,cnt,qd.s,qd.r,xx.wt,qd.wt,prec,maxiter,alpha,skip.iter)
{
    nobs <- dim(r)[1]
    nxi <- dim(r)[2]
    nqd <- dim(qd.r[[1]])[1]
    nx <- length(xx.wt)
    if (!is.null(s)) nnull <- dim(s)[2]
    else nnull <- 0
    nn <- nxi + nnull
    if (is.null(cnt)) cnt <- 0
    ## cv functions
    cv.s <- function(lambda) {
        fit <- .Fortran("cdennewton",
                        cd=as.double(cd), as.integer(nn),
                        as.double(10^(lambda)*q.wk), as.integer(nxi),
                        as.double(t(cbind(r.wk,s))), as.integer(nobs),
                        as.integer(sum(cnt)), as.integer(cnt),
                        as.double(qd.r.wk), as.integer(nqd), as.integer(nx),
                        as.double(xx.wt), as.double(qd.wt),
                        as.double(prec), as.integer(maxiter),
                        as.double(.Machine$double.eps), integer(nn),
                        wk=double(2*(nqd+1)*nx+2*nobs+nn*(2*nn+5)),
                        info=integer(1),PACKAGE="gss")
        if (fit$info==1) stop("gss error in sscden: Newton iteration diverges")
        if (fit$info==2) warning("gss warning in sscden: Newton iteration fails to converge")
        assign("eta",fit$wk[1:(nqd*nx)],inherits=TRUE)
        assign("cd",fit$cd,inherits=TRUE)
        cv <- alpha*fit$wk[nqd*nx+2]-fit$wk[nqd*nx+1]
        alpha.wk <- max(0,log.la0-lambda[1]-5)*(3-alpha) + alpha
        alpha.wk <- min(alpha.wk,3)
        adj <- ifelse (alpha.wk>alpha,(alpha.wk-alpha)*fit$wk[nqd*nx+2],0)
        cv+adj
    }
    cv.m <- function(theta) {
        ind.wk <- theta[1:nq]!=theta.old
        if (sum(ind.wk)==nq) {
            r.wk0 <- 0
            qd.r.wk0 <- array(0,c(nqd,nxi,nx))
            for (i in 1:nq) {
                r.wk0 <- r.wk0 + 10^theta[i]*r[,,i]
                if (length(dim(qd.r[[i]]))==3) qd.r.wk0 <- qd.r.wk0 + 10^theta[i]*qd.r[[i]]
                else qd.r.wk0 <- qd.r.wk0 + as.vector(10^theta[i]*qd.r[[i]])
            }
            assign("r.wk",r.wk0+0,inherits=TRUE)
            assign("qd.r.wk",qd.r.wk0+0,inherits=TRUE)
            assign("theta.old",theta[1:nq]+0,inherits=TRUE)
        }
        else {
            r.wk0 <- r.wk
            qd.r.wk0 <- qd.r.wk
            for (i in (1:nq)[ind.wk]) {
                theta.wk <- (10^(theta[i]-theta.old[i])-1)*10^theta.old[i]
                r.wk0 <- r.wk0 + theta.wk*r[,,i]
                if (length(dim(qd.r[[i]]))==3) qd.r.wk0 <- qd.r.wk0 + theta.wk*qd.r[[i]]
                else qd.r.wk0 <- qd.r.wk0 + as.vector(theta.wk*qd.r[[i]])
            }
        }
        q.wk <- 10^(lambda)*r.wk0[id.basis,]
        qd.r.wk0 <- aperm(qd.r.wk0,c(1,3,2))
        qd.r.wk0 <- array(c(qd.r.wk0,qd.s),c(nqd,nx,nn))
        qd.r.wk0 <- aperm(qd.r.wk0,c(1,3,2))
        fit <- .Fortran("cdennewton",
                        cd=as.double(cd), as.integer(nn),
                        as.double(q.wk), as.integer(nxi),
                        as.double(t(cbind(r.wk0,s))), as.integer(nobs),
                        as.integer(sum(cnt)), as.integer(cnt),
                        as.double(qd.r.wk0), as.integer(nqd), as.integer(nx),
                        as.double(xx.wt), as.double(qd.wt),
                        as.double(prec), as.integer(maxiter),
                        as.double(.Machine$double.eps), integer(nn),
                        wk=double(2*(nqd+1)*nx+2*nobs+nn*(2*nn+5)),
                        info=integer(1),PACKAGE="gss")
        if (fit$info==1) stop("gss error in sscden: Newton iteration diverges")
        if (fit$info==2) warning("gss warning in sscden: Newton iteration fails to converge")
        assign("eta",fit$wk[1:(nqd*nx)],inherits=TRUE)
        assign("cd",fit$cd,inherits=TRUE)
        cv <- alpha*fit$wk[nqd*nx+2]-fit$wk[nqd*nx+1]
        alpha.wk <- max(0,theta[1:nq]-log.th0-5)*(3-alpha) + alpha
        alpha.wk <- min(alpha.wk,3)
        adj <- ifelse (alpha.wk>alpha,(alpha.wk-alpha)*fit$wk[nqd*nx+2],0)
        cv+adj
    }
    cv.m.wk <- function(theta) cv.scale*cv.m(theta)+cv.shift
    ## Initialization
    theta <- -log10(apply(r[id.basis,,,drop=FALSE],3,function(x)sum(diag(x))))
    nq <- length(theta)
    qd.r.wk <- array(0,c(nqd,nxi,nx))
    for (i in 1:nq) {
        if (length(dim(qd.r[[i]]))==3) qd.r.wk <- qd.r.wk + 10^theta[i]*qd.r[[i]]
        else qd.r.wk <- qd.r.wk + as.vector(10^theta[i]*qd.r[[i]])
    }
    if (!nnull) {
        vv.r <- 0
        for (i in 1:nx) {
            mu.r <- apply(qd.r.wk[,,i,drop=FALSE],2,sum)/nqd
            v.r <- apply(qd.r.wk[,,i,drop=FALSE]^2,2,sum)/nqd
            v.r <- v.r - mu.r^2
            vv.r <- vv.r + xx.wt[i]*v.r
        }
        theta.wk <- 0
    }
    else {
        vv.s <- vv.r <- 0
        for (i in 1:nx) {
            mu.s <- apply(qd.s[,i,,drop=FALSE],2,sum)/nqd
            v.s <- apply(qd.s[,i,,drop=FALSE]^2,2,sum)/nqd
            v.s <- v.s - mu.s^2
            mu.r <- apply(qd.r.wk[,,i,drop=FALSE],2,sum)/nqd
            v.r <- apply(qd.r.wk[,,i,drop=FALSE]^2,2,sum)/nqd
            v.r <- v.r - mu.r^2
            vv.s <- vv.s + xx.wt[i]*v.s
            vv.r <- vv.r + xx.wt[i]*v.r
        }
        theta.wk <- log10(sum(vv.s)/nnull/sum(vv.r)*nxi) / 2
    }
    theta <- theta + theta.wk
    qd.r.wk <- aperm(10^theta.wk*qd.r.wk,c(1,3,2))
    qd.r.wk <- array(c(qd.r.wk,qd.s),c(nqd,nx,nn))
    qd.r.wk <- aperm(qd.r.wk,c(1,3,2))
    r.wk <- 0
    for (i in 1:nq) {
        r.wk <- r.wk + 10^theta[i]*r[,,i]
    }
    q.wk <- r.wk[id.basis,]
    log.la0 <- log10(sum(vv.r)/sum(diag(q.wk))) + 2*theta.wk
    ## fixed theta iteration
    eta <- NULL
    cd <- rep(0,nn)
    la <- log.la0
    mn0 <- log.la0-6
    mx0 <- log.la0+6
    repeat {
        mn <- max(la-1,mn0)
        mx <- min(la+1,mx0)
        zz <- nlm0(cv.s,c(mn,mx))
        if ((min(zz$est-mn,mx-zz$est)>=1e-1)||
            (min(zz$est-mn0,mx0-zz$est)<1e-1)) break
        else la <- zz$est
    }
    if (nq==1) {
        jk1 <- cv.s(zz$est)
        c <- cd[1:nxi]
        if (nnull) d <- cd[nxi+(1:nnull)]
        else d <- NULL
        eta <- matrix(eta,nqd,nx)
        for (i in 1:nx) eta[,i] <- eta[,i]/sum(eta[,i])
        return(list(lambda=zz$est,theta=theta,c=c,d=d,cv=jk1,fit=t(eta)))
    }
    ## theta adjustment
    qd.r.wk <- array(0,c(nqd,nxi,nx))
    for (i in 1:nq) {
        theta[i] <- 2*theta[i] + log10(t(cd[1:nxi])%*%r[id.basis,,i]%*%cd[1:nxi])
        if (length(dim(qd.r[[i]]))==3) qd.r.wk <- qd.r.wk + 10^theta[i]*qd.r[[i]]
        else qd.r.wk <- qd.r.wk + as.vector(10^theta[i]*qd.r[[i]])
    }
    if (!nnull) {
        vv.r <- 0
        for (i in 1:nx) {
            mu.r <- apply(qd.r.wk[,,i,drop=FALSE],2,sum)/nqd
            v.r <- apply(qd.r.wk[,,i,drop=FALSE]^2,2,sum)/nqd
            v.r <- v.r - mu.r^2
            vv.r <- vv.r + xx.wt[i]*v.r
        }
        theta.wk <- 0
    }
    else {
        vv.s <- vv.r <- 0
        for (i in 1:nx) {
            mu.s <- apply(qd.s[,i,,drop=FALSE],2,sum)/nqd
            v.s <- apply(qd.s[,i,,drop=FALSE]^2,2,sum)/nqd
            v.s <- v.s - mu.s^2
            mu.r <- apply(qd.r.wk[,,i,drop=FALSE],2,sum)/nqd
            v.r <- apply(qd.r.wk[,,i,drop=FALSE]^2,2,sum)/nqd
            v.r <- v.r - mu.r^2
            vv.s <- vv.s + xx.wt[i]*v.s
            vv.r <- vv.r + xx.wt[i]*v.r
        }
        theta.wk <- log10(sum(vv.s)/nnull/sum(vv.r)*nxi) / 2
    }
    theta <- theta + theta.wk
    qd.r.wk <- aperm(10^theta.wk*qd.r.wk,c(1,3,2))
    qd.r.wk <- array(c(qd.r.wk,qd.s),c(nqd,nx,nn))
    qd.r.wk <- aperm(qd.r.wk,c(1,3,2))
    r.wk <- 0
    for (i in 1:nq) {
        r.wk <- r.wk + 10^theta[i]*r[,,i]
    }
    q.wk <- r.wk[id.basis,]
    log.la0 <- log10(sum(vv.r)/sum(diag(q.wk))) + 2*theta.wk
    log.th0 <- theta-log.la0
    ## fixed theta iteration
    cd <- rep(0,nn)
    la <- log.la0
    mn0 <- log.la0-6
    mx0 <- log.la0+6
    repeat {
        mn <- max(la-1,mn0)
        mx <- min(la+1,mx0)
        zz <- nlm0(cv.s,c(mn,mx))
        if ((min(zz$est-mn,mx-zz$est)>=1e-1)||
            (min(zz$est-mn0,mx0-zz$est)<1e-1)) break
        else la <- zz$est
    }
    lambda <- zz$est
    ## early return
    if (skip.iter) {
        jk1 <- cv.s(zz$est)
        c <- cd[1:nxi]
        if (nnull) d <- cd[nxi+(1:nnull)]
        else d <- NULL
        eta <- matrix(eta,nqd,nx)
        for (i in 1:nx) eta[,i] <- eta[,i]/sum(eta[,i])
        return(list(lambda=zz$est,theta=theta,c=c,d=d,cv=jk1,fit=t(eta)))
    }
    ## theta search
    counter <- 0
    r.wk <- 0
    qd.r.wk <- array(0,c(nqd,nxi,nx))
    for (i in 1:nq) {
        r.wk <- r.wk + 10^theta[i]*r[,,i]
        if (length(dim(qd.r[[i]]))==3) qd.r.wk <- qd.r.wk + 10^theta[i]*qd.r[[i]]
        else qd.r.wk <- qd.r.wk + as.vector(10^theta[i]*qd.r[[i]])
    }
    theta.old <- theta
    tmp <- abs(cv.m(theta))
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
        zz <- nlm(cv.m.wk,theta,stepmax=1,ndigit=7)
        if (zz$code<=3)  break
        theta <- zz$est
        counter <- counter + 1
        if (counter>=5) {
            warning("gss warning in sscden: CV iteration fails to converge")
            break
        }
    }
    ## return
    jk1 <- cv.m(zz$est)
    c <- cd[1:nxi]
    if (nnull) d <- cd[nxi+(1:nnull)]
    else d <- NULL
    eta <- matrix(eta,nqd,nx)
    for (i in 1:nx) eta[,i] <- eta[,i]/sum(eta[,i])
    return(list(lambda=lambda,theta=zz$est,c=c,d=d,cv=jk1,fit=t(eta)))
}
