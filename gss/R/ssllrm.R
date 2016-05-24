## Fit log-linear regression model
ssllrm <- function(formula,response,type=NULL,data=list(),weights,
                   subset,na.action=na.omit,alpha=1,
                   id.basis=NULL,nbasis=NULL,seed=NULL,random=NULL,
                   prec=1e-7,maxiter=30,skip.iter=FALSE)
{
    ## Obtain model frame and model terms
    mf <- match.call()
    mf$response <- mf$type <- mf$alpha <- NULL
    mf$id.basis <- mf$nbasis <- mf$seed <- NULL
    mf$random <- mf$prec <- mf$maxiter <- mf$skip.iter <- NULL
    term.wk <- terms.formula(formula)
    ynames <- as.character(attr(terms(response),"variables"))[-1]
    mf[[1]] <- as.name("model.frame")
    mf <- eval(mf,parent.frame())
    cnt <- model.weights(mf)
    mf$"(weights)" <- NULL
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
            stop("gss error in ssllrm: id.basis out of range")
        nbasis <- length(id.basis)
    }
    ## Check inputs
    mt <- attr(mf,"terms")
    vars <- as.character(attr(mt,"variables"))[-1]
    if(!all(ynames%in%vars)) stop("gss error in ssllrm: response missing in model")
    for (ylab in ynames) {
        if (!is.factor(mf[,ylab])) stop("gss error in ssllrm: response not a factor")
    }
    xnames <- vars[!(vars%in%ynames)]
    if (is.null(xnames)) stop("gss error in ssllrm: missing covariate")
    ## Generate terms    
    term <- mkterm(mf,type)
    term.labels <- labels(mt)
    facs <- attr(mt,"factors")
    ind.wk <- NULL
    for (lab in term.labels)
        ind.wk <- c(ind.wk,any(facs[ynames,lab]))
    term$labels <- term.labels[ind.wk]
    ## Generate quadrature
    qd.pt <- data.frame(levels(mf[,ynames[1]]))
    if (length(ynames)>1) {
        for (ylab in ynames[-1]) {
            wk <- expand.grid(levels(mf[,ylab]),1:dim(qd.pt)[1])
            qd.pt <- data.frame(qd.pt[wk[,2],],wk[,1])
        }
    }
    colnames(qd.pt) <- ynames
    nmesh <- dim(qd.pt)[1]
    x <- mf[,xnames,drop=FALSE]
    ## obtain unique covariate observations
    xx <- mf[,xnames,drop=FALSE]
    if (!is.null(random)) {
        if (class(random)=="formula") random <- mkran(random,data)
        xx <- cbind(xx,random$z)
    }
    xx <- apply(xx,1,function(x)paste(x,collapse="\r"))
    x.dup.ind <- duplicated(xx)    
    if (!is.null(cnt)) xx <- rep(xx,cnt)
    xx.wt <- as.vector(table(xx)[unique(xx)])
    xx.wt <- xx.wt/sum(xx.wt)
    nx <- length(xx.wt)
    ## Generate Random
    if (!is.null(random)) {
        ## z and qd.z
        z <- qd.z <- nlvl <- NULL
        for (ylab in ynames) {
            y.wk <- mf[,ylab]
            pt.wk <- qd.pt[,ylab]
            lvl.wk <- levels(y.wk)
            nlvl.wk <- length(lvl.wk)
            nlvl <- c(nlvl,nlvl.wk)
            z.aux <- diag(1,nlvl.wk-1)
            z.aux <- rbind(z.aux,rep(-1,nlvl.wk-1))
            rownames(z.aux) <- lvl.wk
            for (i in 1:(nlvl.wk-1)) {
                z <- cbind(z,z.aux[y.wk,i]*random$z)
                for (j in 1:nmesh) {
                    qd.z <- cbind(qd.z,z.aux[pt.wk[j],i]*random$z[!x.dup.ind,])
                }
            }
        }
        nz <- dim(random$z)[2]
        nZ <- sum(nlvl-1)*nz
        qd.z <- aperm(array(qd.z,c(nx,nz,nmesh,nZ/nz)),c(3,1,2,4))
        qd.z <- array(qd.z,c(nmesh,nx,nZ))
        ## Sigma
        env <- list(sigma=random$sigma,nzeta=length(random$init),nz=nz,nlvl=nlvl)
        fun <- function(zeta,env) {
            ny <- length(env$nlvl)
            nze <- env$nzeta
            sigma <- env$sigma
            dm <- cumsum(env$nlvl-1)*env$nz
            zz <- matrix(0,dm[ny],dm[ny])
            dm <- c(0,dm)
            for (i in 1:ny) {
                nlvl.wk <- nlvl[i]
                wk <- kronecker(diag(1,nlvl.wk-1)+1,
                                sigma$fun(zeta[nze*(i-1)+(1:nze)],sigma$env))
                zz[(dm[i]+1):dm[i+1],(dm[i]+1):dm[i+1]] <- wk
            }
            zz
        }
        Sigma <- list(fun=fun,env=env)
        ## init
        init <- rep(random$init,length(nlvl))
        ## assemble
        Random <- list(z=z,qd.z=qd.z,sigma=Sigma,init=init)
    }
    else Random <- NULL
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
            stop("gss error in ssllrm: unpenalized terms are linearly dependent")
    }
    ## Fit the model
    z <- mspllrm(s,r,id.basis,cnt,qd.s,qd.r,xx.wt,Random,prec,maxiter,alpha,skip.iter)
    ## Brief description of model terms
    desc <- NULL
    for (label in term$labels)
        desc <- rbind(desc,as.numeric(c(term[[label]][c("nphi","nrk")])))
    desc <- rbind(desc,apply(desc,2,sum))
    rownames(desc) <- c(term$labels,"total")
    colnames(desc) <- c("Unpenalized","Penalized")
    ## Return the results
    obj <- c(list(call=match.call(),mf=mf,cnt=cnt,terms=term,desc=desc,
                  qd.pt=qd.pt,xx.wt=xx.wt,x.dup.ind=x.dup.ind,alpha=alpha,
                  ynames=ynames,xnames=xnames,id.basis=id.basis,
                  random=random,Random=Random,skip.iter=skip.iter),z)
    if (is.null(cnt)) obj$se.aux$v <- sqrt(nobs)*obj$se.aux$v
    else obj$se.aux$v <- sqrt(sum(cnt))*obj$se.aux$v
    class(obj) <- c("ssllrm")
    obj
}

## Fit (multiple smoothing parameter) log-linear regression model
mspllrm <- function(s,r,id.basis,cnt,qd.s,qd.r,xx.wt,random,prec,maxiter,alpha,skip.iter)
{
    nobs <- dim(r)[1]
    nxi <- dim(r)[2]
    nqd <- dim(qd.r[[1]])[1]
    nx <- length(xx.wt)
    if (!is.null(s)) nnull <- dim(s)[2]
    else nnull <- 0
    if (!is.null(random)) nz <- ncol(as.matrix(random$z))
    else nz <- 0
    nxiz <- nxi + nz
    nn <- nxiz + nnull
    if (is.null(cnt)) cnt <- 0
    ## cv functions
    cv.s <- function(lambda) {
        if (is.null(random)) q.wk0 <- 10^(lambda)*q.wk
        else {
            q.wk0 <- matrix(0,nxiz,nxiz)
            q.wk0[1:nxi,1:nxi] <- 10^(lambda[1])*q.wk
            q.wk0[(nxi+1):nxiz,(nxi+1):nxiz] <-
                10^(2*ran.scal)*random$sigma$fun(lambda[-1],random$sigma$env)
        }
        fit <- .Fortran("llrmnewton",
                        cd=as.double(cd), as.integer(nn),
                        as.double(q.wk0), as.integer(nxiz),
                        as.double(t(cbind(r.wk,s))), as.integer(nobs),
                        as.integer(sum(cnt)), as.integer(cnt),
                        as.double(qd.r.wk), as.integer(nqd), as.integer(nx),
                        as.double(xx.wt),
                        as.double(prec), as.integer(maxiter),
                        as.double(.Machine$double.eps), integer(nn),
                        wk=double(2*(nqd+1)*nx+2*nobs+nn*(2*nn+5)),
                        info=integer(1),PACKAGE="gss")
        if (fit$info==1) stop("gss error in ssllrm: Newton iteration diverges")
        if (fit$info==2) warning("gss warning in ssllrm: Newton iteration fails to converge")
        assign("eta",fit$wk[1:(nqd*nx)],inherits=TRUE)
        assign("cd",fit$cd,inherits=TRUE)
        cv <- alpha*fit$wk[nqd*nx+2]-fit$wk[nqd*nx+1]
        alpha.wk <- max(0,log.la0-lambda[1]-5)*(3-alpha) + alpha
        alpha.wk <- min(alpha.wk,3)
        adj <- ifelse (alpha.wk>alpha,(alpha.wk-alpha)*fit$wk[nqd*nx+2],0)
        cv+adj
    }
    cv.s.wk <- function(lambda) cv.scale*cv.s(lambda)+cv.shift
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
        q.wk <- r.wk0[id.basis,]
        if (is.null(random)) q.wk0 <- 10^(lambda)*q.wk
        else {
            r.wk0 <- cbind(r.wk0,10^ran.scal*random$z)
            q.wk0 <- matrix(0,nxiz,nxiz)
            q.wk0[1:nxi,1:nxi] <- 10^(lambda)*q.wk
            q.wk0[(nxi+1):nxiz,(nxi+1):nxiz] <-
              10^(2*ran.scal)*random$sigma$fun(theta[-(1:nq)],random$sigma$env)
        }
        qd.r.wk0 <- aperm(qd.r.wk0,c(1,3,2))
        if (!is.null(random)) {
            qd.r.wk0 <- array(c(qd.r.wk0,10^ran.scal*random$qd.z),c(nqd,nx,nxiz))
        }
        qd.r.wk0 <- array(c(qd.r.wk0,qd.s),c(nqd,nx,nn))
        qd.r.wk0 <- aperm(qd.r.wk0,c(1,3,2))
        fit <- .Fortran("llrmnewton",
                        cd=as.double(cd), as.integer(nn),
                        as.double(q.wk0), as.integer(nxiz),
                        as.double(t(cbind(r.wk0,s))), as.integer(nobs),
                        as.integer(sum(cnt)), as.integer(cnt),
                        as.double(qd.r.wk0), as.integer(nqd), as.integer(nx),
                        as.double(xx.wt),
                        as.double(prec), as.integer(maxiter),
                        as.double(.Machine$double.eps), integer(nn),
                        wk=double(2*(nqd+1)*nx+2*nobs+nn*(2*nn+5)),
                        info=integer(1),PACKAGE="gss")
        if (fit$info==1) stop("gss error in ssllrm: Newton iteration diverges")
        if (fit$info==2) warning("gss warning in ssllrm: Newton iteration fails to converge")
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
    if (!is.null(random)) {
        vv.z <- 0
        for (i in 1:nx) {
            mu.z <- apply(random$qd.z[,i,,drop=FALSE],2,sum)/nqd
            v.z <- apply(random$qd.z[,i,,drop=FALSE]^2,2,sum)/nqd
            v.z <- v.z - mu.z^2
            vv.z <- vv.z + xx.wt[i]*v.z
        }
        ran.scal <- theta.wk - log10(sum(vv.z)/nz/sum(vv.r)*nxi) / 2
    }
    else ran.scal <- NULL
    theta <- theta + theta.wk
    qd.r.wk <- aperm(10^theta.wk*qd.r.wk,c(1,3,2))
    if (!is.null(random)) {
        qd.r.wk <- array(c(qd.r.wk,10^ran.scal*random$qd.z),c(nqd,nx,nxiz))
    }
    qd.r.wk <- array(c(qd.r.wk,qd.s),c(nqd,nx,nn))
    qd.r.wk <- aperm(qd.r.wk,c(1,3,2))
    r.wk <- 0
    for (i in 1:nq) {
        r.wk <- r.wk + 10^theta[i]*r[,,i]
    }
    q.wk <- r.wk[id.basis,]
    if (!is.null(random)) r.wk <- cbind(r.wk,10^ran.scal*random$z)
    log.la0 <- log10(sum(vv.r)/sum(diag(q.wk))) + 2*theta.wk
    ## fixed theta iteration
    eta <- NULL
    cd <- rep(0,nn)
    if (is.null(random)) la <- log.la0
    else la <- c(log.la0,random$init)
    if (length(la)-1) {
        counter <- 0
        ## scale and shift cv
        tmp <- abs(cv.s(la))
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
            zz <- nlm(cv.s.wk,la,stepmax=1,ndigit=7)
            if (zz$code<=3) break
            la <- zz$est
            counter <- counter + 1
            if (counter>=5) {
                warning("gss warning in ssllrm: iteration for model selection fails to converge")
                break
            }
        }
        cv <- (zz$min-cv.shift)/cv.scale
    }
    else {
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
        cv <- zz$min
    }
    if (nq==1) {
        if (is.null(random)) {
            lambda <- zz$est
            zeta <- NULL
        }
        else {
            lambda <- zz$est[1]
            zeta <- zz$est[-1]
        }
        if (is.null(random)) q.wk0 <- 10^(lambda)*q.wk
        else {
            q.wk0 <- matrix(0,nxiz,nxiz)
            q.wk0[1:nxi,1:nxi] <- 10^(lambda)*q.wk
            q.wk0[(nxi+1):nxiz,(nxi+1):nxiz] <-
                10^(2*ran.scal)*random$sigma$fun(zeta,random$sigma$env)
        }
        se.aux <- .Fortran("llrmaux",
                           as.double(cd), as.integer(nn),
                           as.double(q.wk0), as.integer(nxiz),
                           as.double(qd.r.wk), as.integer(nqd),
                           as.integer(nx), as.double(xx.wt),
                           as.double(.Machine$double.eps), double(nqd*nx),
                           double(nx), double(nn),
                           v=double(nn*nn), double(nn*nn),
                           jpvt=integer(nn), PACKAGE="gss")[c("v","jpvt")]
        c <- cd[1:nxi]
        if (nz) b <- 10^ran.scal*cd[nxi+(1:nz)]
        else b <- NULL
        if (nnull) d <- cd[nxiz+(1:nnull)]
        else d <- NULL
        eta <- matrix(eta,nqd,nx)
        for (i in 1:nx) eta[,i] <- eta[,i]/sum(eta[,i])
        return(list(lambda=lambda,zeta=zeta,theta=theta,ran.scal=ran.scal,
                    c=c,b=b,d=d,cv=cv,fit=t(eta),se.aux=se.aux))
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
    if (!is.null(random)) {
        vv.z <- 0
        for (i in 1:nx) {
            mu.z <- apply(random$qd.z[,i,,drop=FALSE],2,sum)/nqd
            v.z <- apply(random$qd.z[,i,,drop=FALSE]^2,2,sum)/nqd
            v.z <- v.z - mu.z^2
            vv.z <- vv.z + xx.wt[i]*v.z
        }
        ran.scal <- theta.wk - log10(sum(vv.z)/nz/sum(vv.r)*nxi) / 2
    }
    theta <- theta + theta.wk
    qd.r.wk <- aperm(10^theta.wk*qd.r.wk,c(1,3,2))
    if (!is.null(random)) {
        qd.r.wk <- array(c(qd.r.wk,10^ran.scal*random$qd.z),c(nqd,nx,nxiz))
    }
    qd.r.wk <- array(c(qd.r.wk,qd.s),c(nqd,nx,nn))
    qd.r.wk <- aperm(qd.r.wk,c(1,3,2))
    r.wk <- 0
    for (i in 1:nq) {
        r.wk <- r.wk + 10^theta[i]*r[,,i]
    }
    q.wk <- r.wk[id.basis,]
    if (!is.null(random)) r.wk <- cbind(r.wk,10^ran.scal*random$z)
    log.la0 <- log10(sum(vv.r)/sum(diag(q.wk))) + 2*theta.wk
    log.th0 <- theta-log.la0
    ## fixed theta iteration
    cd <- rep(0,nn)
    if (is.null(random)) la <- log.la0
    else la <- c(log.la0,random$init)
    if (length(la)-1) {
        counter <- 0
        ## scale and shift cv
        tmp <- abs(cv.s(la))
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
            zz <- nlm(cv.s.wk,la,stepmax=1,ndigit=7)
            if (zz$code<=3) break
            la <- zz$est
            counter <- counter + 1
            if (counter>=5) {
                warning("gss warning in ssllrm: iteration for model selection fails to converge")
                break
            }
        }
        cv <- (zz$min-cv.shift)/cv.scale
    }
    else {
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
        cv <- zz$min
    }
    if (is.null(random)) {
        lambda <- zz$est
        zeta <- NULL
    }
    else {
        lambda <- zz$est[1]
        zeta <- zz$est[-1]
    }
    ## early return
    if (skip.iter) {
        if (is.null(random)) q.wk0 <- 10^(lambda)*q.wk
        else {
            q.wk0 <- matrix(0,nxiz,nxiz)
            q.wk0[1:nxi,1:nxi] <- 10^(lambda)*q.wk
            q.wk0[(nxi+1):nxiz,(nxi+1):nxiz] <-
                10^(2*ran.scal)*random$sigma$fun(zeta,random$sigma$env)
        }
        se.aux <- .Fortran("llrmaux",
                           as.double(cd), as.integer(nn),
                           as.double(q.wk0), as.integer(nxiz),
                           as.double(qd.r.wk), as.integer(nqd),
                           as.integer(nx), as.double(xx.wt),
                           as.double(.Machine$double.eps), double(nqd*nx),
                           double(nx), double(nn),
                           v=double(nn*nn), double(nn*nn),
                           jpvt=integer(nn), PACKAGE="gss")[c("v","jpvt")]
        c <- cd[1:nxi]
        if (nz) b <- 10^ran.scal*cd[nxi+(1:nz)]
        else b <- NULL
        if (nnull) d <- cd[nxiz+(1:nnull)]
        else d <- NULL
        eta <- matrix(eta,nqd,nx)
        for (i in 1:nx) eta[,i] <- eta[,i]/sum(eta[,i])
        return(list(lambda=lambda,zeta=zeta,theta=theta,ran.scal=ran.scal,
                    c=c,b=b,d=d,cv=cv,fit=t(eta),se.aux=se.aux))
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
    if (!is.null(random)) theta <- c(theta,zeta)
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
            warning("gss warning in ssllrm: CV iteration fails to converge")
            break
        }
    }
    cv <- (zz$min-cv.shift)/cv.scale
    if (is.null(random)) {
        theta <- zz$est
        zeta <- NULL
    }
    else {
        theta <- zz$est[1:nq]
        zeta <- zz$est[-(1:nq)]
    }
    ## return
    q.wk <- 0
    qd.r.wk <- array(0,c(nqd,nxi,nx))
    for (i in 1:nq) {
        q.wk <- q.wk + 10^theta[i]*r[id.basis,,i]
        if (length(dim(qd.r[[i]]))==3) qd.r.wk <- qd.r.wk + 10^theta[i]*qd.r[[i]]
        else qd.r.wk <- qd.r.wk + as.vector(10^theta[i]*qd.r[[i]])
    }
    if (is.null(random)) q.wk0 <- 10^(lambda)*q.wk
    else {
        q.wk0 <- matrix(0,nxiz,nxiz)
        q.wk0[1:nxi,1:nxi] <- 10^(lambda)*q.wk
        q.wk0[(nxi+1):nxiz,(nxi+1):nxiz] <-
          10^(2*ran.scal)*random$sigma$fun(zeta,random$sigma$env)
    }
    qd.r.wk <- aperm(qd.r.wk,c(1,3,2))
    if (!is.null(random)) {
        qd.r.wk <- array(c(qd.r.wk,10^ran.scal*random$qd.z),c(nqd,nx,nxiz))
    }
    qd.r.wk <- array(c(qd.r.wk,qd.s),c(nqd,nx,nn))
    qd.r.wk <- aperm(qd.r.wk,c(1,3,2))
    se.aux <- .Fortran("llrmaux",
                       as.double(cd), as.integer(nn),
                       as.double(q.wk0), as.integer(nxiz),
                       as.double(qd.r.wk), as.integer(nqd),
                       as.integer(nx), as.double(xx.wt),
                       as.double(.Machine$double.eps), double(nqd*nx),
                       double(nx), double(nn),
                       v=double(nn*nn), double(nn*nn),
                       jpvt=integer(nn), PACKAGE="gss")[c("v","jpvt")]
    c <- cd[1:nxi]
    if (nz) b <- 10^ran.scal*cd[nxi+(1:nz)]
    else b <- NULL
    if (nnull) d <- cd[nxiz+(1:nnull)]
    else d <- NULL
    eta <- matrix(eta,nqd,nx)
    for (i in 1:nx) eta[,i] <- eta[,i]/sum(eta[,i])
    return(list(lambda=lambda,zeta=zeta,theta=theta,ran.scal=ran.scal,
                c=c,b=b,d=d,cv=cv,fit=t(eta),se.aux=se.aux))
}
