## Fit hazard model
sscox <- function(formula,type=NULL,data=list(),weights=NULL,subset,
                  na.action=na.omit,partial=NULL,alpha=1.4,
                  id.basis=NULL,nbasis=NULL,seed=NULL,random=NULL,
                  prec=1e-7,maxiter=30,skip.iter=FALSE)
{
    ## Local functions handling formula
    Surv <- function(time,status,start=0) {
        if (!is.numeric(time)|!is.vector(time))
            stop("gss error in sscox: time should be a numerical vector")
        if ((nobs <- length(time))-length(status))
            stop("gss error in sscox: time and status mismatch in size")
        if ((length(start)-nobs)&(length(start)-1))
            stop("gss error in sscox: time and start mismatch in size")
        if (any(start>time))
            stop("gss error in sscox: start after follow-up time")
        if (min(start)<0)
            warning("gss warning in sscox: start before time 0")
        time <- cbind(start,time)
        list(start=time[,1],end=time[,2],status=as.logical(status))
    }
    ## Obtain model frame and model terms
    mf <- match.call()
    mf$type <- mf$alpha <- mf$random <- mf$partial <- NULL
    mf$id.basis <- mf$nbasis <- mf$seed <- NULL
    mf$prec <- mf$maxiter <- mf$skip.iter <- NULL
    term.wk <- terms.formula(formula)
    ## response
    resp <- attr(term.wk,"variable")[[2]]
    ind.wk <- length(strsplit(deparse(resp),'')[[1]])
    if ((substr(deparse(resp),1,5)!='Surv(')
        |(substr(deparse(resp),ind.wk,ind.wk)!=')'))
        stop("gss error in sscox: response should be Surv(...)")
    yy <- with(data,eval(resp))
    ## model frame
    term.labels <- attr(term.wk,"term.labels")
    mf[[1]] <- as.name("model.frame")
    mf[[2]] <- eval(parse(text=paste("~",paste(term.labels,collapse="+"),"-1")))
    mf <- eval(mf,parent.frame())
    ## trim yy if subset is used
    nobs <- nrow(mf)
    if (nobs<length(yy$end)) {
        yy$start <- yy$start[subset]
        yy$end <- yy$end[subset]
        yy$status <- yy$status[subset]
    }
    ## Generate sub-basis
    cnt <- model.weights(mf)
    if (!is.null(cnt)) mf["(weights)"] <- NULL
    if (is.null(id.basis)) {
        if (is.null(nbasis)) nbasis <- max(30,ceiling(10*nobs^(2/9)))
        if (nbasis>sum(yy$status)) nbasis <- sum(yy$status)
        if (!is.null(seed)) set.seed(seed)
        id.basis <- sample((1:nobs)[yy$status],nbasis,prob=cnt[yy$status])
    }
    else {
        if (!all(id.basis%in%(1:nobs)[yy$status]))
            stop("gss error in sscox: id.basis not all at failure cases")
        nbasis <- length(id.basis)
    }
    id.wk <- NULL
    nT <- sum(yy$status)
    for (i in 1:nbasis) {
        id.wk <- c(id.wk,(1:nT)[(1:nobs)[yy$status]%in%id.basis[i]])
    }
    ## Generate terms    
    term <- mkterm(mf,type)
    term$labels <- term$labels[term$labels!="1"]
    ## Generate random
    if (!is.null(random)) {
        if (class(random)=="formula") random <- mkran(random,data)
        random$qd.z <- random$z
        random$z <- random$z[yy$status,]
    }
    ## Generate s and r
    s <- qd.s <- r <- qd.r <- NULL
    nq <- 0
    for (label in term$labels) {
        x.basis <- mf[id.basis,term[[label]]$vlist]
        qd.x <- mf[,term[[label]]$vlist]
        nphi <- term[[label]]$nphi
        nrk <- term[[label]]$nrk
        if (nphi) {
            phi <- term[[label]]$phi
            for (i in 1:nphi) {
                s.wk <- phi$fun(qd.x,nu=i,env=phi$env)
                s <- cbind(s,s.wk[yy$status])
                qd.s <- cbind(qd.s,s.wk)
            }
        }
        if (nrk) {
            rk <- term[[label]]$rk
            for (i in 1:nrk) {
                nq <- nq+1
                r.wk <- rk$fun(qd.x,x.basis,nu=i,env=rk$env,out=TRUE)
                r <- array(c(r,r.wk[yy$status,]),c(nT,nbasis,nq))
                qd.r <- array(c(qd.r,r.wk),c(nobs,nbasis,nq))
            }
        }
    }
    ## Add the partial term
    if (!is.null(partial)) {
        mf.p <- model.frame(partial,data)
        for (lab in colnames(mf.p)) mf[,lab] <- mf.p[,lab]
        mt.p <- attr(mf.p,"terms")
        lab.p <- labels(mt.p)
        matx.p <- model.matrix(mt.p,data)[,-1,drop=FALSE]
        if (dim(matx.p)[1]!=dim(mf)[1])
            stop("gss error in sscox: partial data are of wrong size")
        matx.p <- scale(matx.p)
        center.p <- attr(matx.p,"scaled:center")
        scale.p <- attr(matx.p,"scaled:scale")
        s <- cbind(s,matx.p[yy$status,])
        qd.s <- cbind(qd.s,matx.p)
        part <- list(mt=mt.p,center=center.p,scale=scale.p)
    }
    else part <- lab.p <- NULL
    ## Check s rank
    if (!is.null(s)) {
        nnull <- dim(s)[2]
        if (qr(s)$rank<nnull)
            stop("gss error in sscox: unpenalized terms are linearly dependent")
    }
    ## Generate quadrature and biasing weights
    if (is.null(cnt)) {
        qd.wt <- rep(1,dim(mf)[1])
        cntt <- NULL
        b.wt <- rep(1/nT,nT)
    }
    else {
        qd.wt <- cnt
        cntt <- cnt[yy$status]
        b.wt <- cntt/sum(cntt)
    }
    tt <- yy$end[yy$status]
    t.wt <- (outer(yy$end,tt,">=")&outer(yy$start,tt,"<="))/1
    bias0 <- list(nt=nT,wt=b.wt,qd.wt=t.wt)
    ## Fit the model
    if (nq==1) {
        r <- r[,,1]
        qd.r <- qd.r[,,1]
        z <- sspcox(s,r,r[id.wk,],cntt,qd.s,qd.r,qd.wt,prec,maxiter,alpha,random,bias0)
    }
    else z <- mspcox(s,r,id.wk,cntt,qd.s,qd.r,qd.wt,prec,maxiter,alpha,random,bias0,skip.iter)
    ## Brief description of model terms
    desc <- NULL
    for (label in term$labels)
        desc <- rbind(desc,as.numeric(c(term[[label]][c("nphi","nrk")])))
    if (!is.null(partial)) {
        desc <- rbind(desc,matrix(c(1,0),length(lab.p),2,byrow=TRUE))
    }
    desc <- rbind(desc,apply(desc,2,sum))
    if (is.null(partial)) rownames(desc) <- c(term$labels,"total")
    else rownames(desc) <- c(term$labels,lab.p,"total")
    colnames(desc) <- c("Unpenalized","Penalized")
    ## Return the results
    obj <- c(list(call=match.call(),mf=mf,cnt=cnt,terms=term,desc=desc,
                  alpha=alpha,id.basis=id.basis,partial=part,lab.p=lab.p,
                  random=random,bias=bias0,skip.iter=skip.iter),z)
    Nobs <- ifelse(is.null(cnt),nT,sum(cntt))
    obj$se.aux$v <- sqrt(Nobs)*obj$se.aux$v
    class(obj) <- c("sscox")
    obj
}

## Fit single smoothing parameter density
sspcox <- function(s,r,q,cnt,qd.s,qd.r,qd.wt,prec,maxiter,alpha,random,bias)
{
    nobs <- dim(r)[1]
    nxi <- dim(r)[2]
    nqd <- length(qd.wt)
    if (!is.null(s)) nnull <- dim(s)[2]
    else nnull <- 0
    if (!is.null(random)) nz <- ncol(as.matrix(random$z))
    else nz <- 0
    nxiz <- nxi + nz
    nn <- nxiz + nnull
    if (is.null(cnt)) cnt <- 0
    ## cv function
    cv <- function(lambda) {
        if (is.null(random)) q.wk0 <- 10^(lambda+theta)*q
        else {
            q.wk0 <- matrix(0,nxiz,nxiz)
            q.wk0[1:nxi,1:nxi] <- 10^(lambda[1]+theta)*q
            q.wk0[(nxi+1):nxiz,(nxi+1):nxiz] <-
                10^(2*ran.scal)*random$sigma$fun(lambda[-1],random$sigma$env)
        }
        fit <- .Fortran("dnewton",
                        cd=as.double(cd), as.integer(nn),
                        as.double(q.wk0), as.integer(nxiz),
                        as.double(t(cbind(r.wk,s))), as.integer(nobs),
                        as.integer(sum(cnt)), as.integer(cnt),
                        as.double(cbind(qd.r.wk,qd.s)), as.integer(nqd),
                        as.integer(bias$nt), as.double(bias$wt),
                        as.double(t(qd.wt*bias$qd.wt)),
                        as.double(prec), as.integer(maxiter),
                        as.double(.Machine$double.eps), integer(nn),
                        wk=double(2*((nqd+1)*bias$nt+nobs)+nn*(2*nn+4)+max(nn,3)),
                        info=integer(1),PACKAGE="gss")
        if (fit$info==1) stop("gss error in sscox: Newton iteration diverges")
        if (fit$info==2) warning("gss warning in sscox: Newton iteration fails to converge")
        assign("cd",fit$cd,inherits=TRUE)
        cv <- alpha*fit$wk[2]-fit$wk[1]
        alpha.wk <- max(0,log.la0-lambda-5)*(3-alpha) + alpha
        alpha.wk <- min(alpha.wk,3)
        adj <- ifelse (alpha.wk>alpha,(alpha.wk-alpha)*fit$wk[2],0)
        cv+adj
    }
    cv.wk <- function(lambda) cv.scale*cv(lambda)+cv.shift
    ## initialization
    if (!nnull) {
        vv.r <- 0
        for (i in 1:bias$nt) {
            wt.wk <- qd.wt*bias$qd.wt[,i]
            mu.r <- apply(wt.wk*qd.r,2,sum)/sum(wt.wk)
            v.r <- apply(wt.wk*qd.r^2,2,sum)/sum(wt.wk)
            v.r <- v.r - mu.r^2
            vv.r <- vv.r + bias$wt[i]*v.r
        }
        theta <- 0
    }
    else {
        vv.s <- vv.r <- 0
        for (i in 1:bias$nt) {
            wt.wk <- qd.wt*bias$qd.wt[,i]
            mu.s <- apply(wt.wk*qd.s,2,sum)/sum(wt.wk)
            v.s <- apply(wt.wk*qd.s^2,2,sum)/sum(wt.wk)
            v.s <- v.s - mu.s^2
            mu.r <- apply(wt.wk*qd.r,2,sum)/sum(wt.wk)
            v.r <- apply(wt.wk*qd.r^2,2,sum)/sum(wt.wk)
            v.r <- v.r - mu.r^2
            vv.s <- vv.s + bias$wt[i]*v.s
            vv.r <- vv.r + bias$wt[i]*v.r
        }
        theta <- log10(sum(vv.s)/nnull/sum(vv.r)*nxi) / 2
    }
    log.la0 <- log10(sum(vv.r)/sum(diag(q))) + theta
    if (!is.null(random)) {
        mu.z <- apply(qd.wt*random$qd.z,2,sum)
        v.z <- apply(qd.wt*random$qd.z^2,2,sum)
        ran.scal <- theta - log10(sum(v.z-mu.z^2)/nz/sum(v.r-mu.r^2)*nxi) / 2
        r.wk <- cbind(10^theta*r,10^ran.scal*random$z)
        qd.r.wk <- cbind(10^theta*qd.r,10^ran.scal*random$qd.z)
    }
    else {
        ran.scal <- NULL
        r.wk <- 10^theta*r
        qd.r.wk <- 10^theta*qd.r
    }
    ## lambda search
    cd <- rep(0,nn)
    if (is.null(random)) la <- log.la0
    else la <- c(log.la0,random$init)
    if (length(la)-1) {
        counter <- 0
        ## scale and shift cv
        tmp <- abs(cv(la))
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
            zz <- nlm(cv.wk,la,stepmax=1,ndigit=7)
            if (zz$code<=3) break
            la <- zz$est
            counter <- counter + 1
            if (counter>=5) {
                warning("gss warning in sscox: iteration for model selection fails to converge")
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
            zz <- nlm0(cv,c(mn,mx))
            if ((min(zz$est-mn,mx-zz$est)>=1e-1)||
                (min(zz$est-mn0,mx0-zz$est)<1e-1)) break
            else la <- zz$est
        }
        cv <- zz$min
    }
    ## return
    if (is.null(random)) {
        lambda <- zz$est
        zeta <- NULL
    }
    else {
        lambda <- zz$est[1]
        zeta <- zz$est[-1]
    }
    if (is.null(random)) {
        q.wk0 <- 10^(lambda+theta)*q
        qd.r.wk <- 10^theta*qd.r
    }
    else {
        q.wk0 <- matrix(0,nxiz,nxiz)
        q.wk0[1:nxi,1:nxi] <- 10^(lambda+theta)*q
        q.wk0[(nxi+1):nxiz,(nxi+1):nxiz] <-
          10^(2*ran.scal)*random$sigma$fun(zeta,random$sigma$env)
        qd.r.wk <- cbind(10^theta*qd.r,10^ran.scal*random$qd.z)
    }
    se.aux <- .Fortran("coxaux",
                       as.double(cd), as.integer(nn),
                       as.double(q.wk0), as.integer(nxiz),
                       as.double(cbind(qd.r.wk,qd.s)), as.integer(nqd),
                       as.integer(bias$nt), as.double(bias$wt),
                       as.double(.Machine$double.eps),
                       as.double(qd.wt*bias$qd.wt),
                       double(nqd*bias$nt), double(bias$nt),
                       double(nn), v=double(nn*nn), double(nn*nn),
                       jpvt=integer(nn), PACKAGE="gss")[c("v","jpvt")]
    c <- cd[1:nxi]
    if (nz) b <- 10^ran.scal*cd[nxi+(1:nz)]
    else b <- NULL
    if (nnull) d <- cd[nxiz+(1:nnull)]
    else d <- NULL
    return(list(lambda=lambda,zeta=zeta,theta=theta,ran.scal=ran.scal,
                c=c,b=b,d=d,cv=cv,se.aux=se.aux))
}

## Fit multiple smoothing parameter density
mspcox <- function(s,r,id.basis,cnt,qd.s,qd.r,qd.wt,prec,maxiter,alpha,random,bias,skip.iter)
{
    nobs <- dim(r)[1]
    nxi <- dim(r)[2]
    nq <- dim(r)[3]
    nqd <- length(qd.wt)
    if (!is.null(s)) nnull <- dim(s)[2]
    else nnull <- 0
    if (!is.null(random)) nz <- ncol(as.matrix(random$z))
    else nz <- 0
    nxiz <- nxi + nz
    nn <- nxiz + nnull
    if (is.null(cnt)) cnt <- 0
    ## cv function
    cv <- function(theta) {
        ind.wk <- theta[1:nq]!=theta.old
        if (sum(ind.wk)==nq) {
            r.wk0 <- qd.r.wk0 <- 0
            for (i in 1:nq) {
                r.wk0 <- r.wk0 + 10^theta[i]*r[,,i]
                qd.r.wk0 <- qd.r.wk0 + 10^theta[i]*qd.r[,,i]
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
                qd.r.wk0 <- qd.r.wk0 + theta.wk*qd.r[,,i]
            }
        }
        q.wk <- r.wk0[id.basis,]
        if (is.null(random)) q.wk0 <- 10^(lambda)*q.wk
        else {
            r.wk0 <- cbind(r.wk0,10^ran.scal*random$z)
            qd.r.wk0 <- cbind(qd.r.wk0,10^ran.scal*random$qd.z)
            q.wk0 <- matrix(0,nxiz,nxiz)
            q.wk0[1:nxi,1:nxi] <- 10^(lambda)*q.wk
            q.wk0[(nxi+1):nxiz,(nxi+1):nxiz] <-
              10^(2*ran.scal)*random$sigma$fun(theta[-(1:nq)],random$sigma$env)
        }
        fit <- .Fortran("dnewton",
                        cd=as.double(cd), as.integer(nn),
                        as.double(q.wk0), as.integer(nxiz),
                        as.double(t(cbind(r.wk0,s))), as.integer(nobs),
                        as.integer(sum(cnt)), as.integer(cnt),
                        as.double(cbind(qd.r.wk0,qd.s)), as.integer(nqd),
                        as.integer(bias$nt), as.double(bias$wt),
                        as.double(t(qd.wt*bias$qd.wt)),
                        as.double(prec), as.integer(maxiter),
                        as.double(.Machine$double.eps), integer(nn),
                        wk=double(2*((nqd+1)*bias$nt+nobs)+nn*(2*nn+4)+max(nn,3)),
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
    theta <- -log10(apply(r[id.basis,,],3,function(x)sum(diag(x))))
    r.wk <- qd.r.wk <- 0
    for (i in 1:nq) {
        r.wk <- r.wk + 10^theta[i]*r[,,i]
        qd.r.wk <- qd.r.wk + 10^theta[i]*qd.r[,,i]
    }
    ## theta adjustment
    z <- sspcox(s,r.wk,r.wk[id.basis,],cnt,qd.s,qd.r.wk,qd.wt,prec,maxiter,alpha,random,bias)
    theta <- theta + z$theta
    r.wk <- qd.r.wk <- 0
    for (i in 1:nq) {
        theta[i] <- 2*theta[i] + log10(t(z$c)%*%r[id.basis,,i]%*%z$c)
        r.wk <- r.wk + 10^theta[i]*r[,,i]
        qd.r.wk <- qd.r.wk + 10^theta[i]*qd.r[,,i]
    }
    mu <- apply(qd.wt*qd.r.wk,2,sum)/sum(qd.wt)
    v <- apply(qd.wt*qd.r.wk^2,2,sum)/sum(qd.wt)
    log.la0 <- log10(sum(v-mu^2)/sum(diag(r.wk[id.basis,])))
    log.th0 <- theta-log.la0
    ## lambda search
    z <- sspcox(s,r.wk,r.wk[id.basis,],cnt,qd.s,qd.r.wk,qd.wt,prec,maxiter,alpha,random,bias)
    ## early return
    if (skip.iter) {
        z$theta <- theta
        return(z)
    }
    ## theta search
    lambda <- z$lambda
    log.th0 <- log.th0 + z$lambda
    theta <- theta + z$theta
    ran.scal <- z$ran.scal
    cd <- c(z$c,z$b,z$d)
    counter <- 0
    r.wk <- qd.r.wk <- 0
    for (i in 1:nq) {
        r.wk <- r.wk + 10^theta[i]*r[,,i]
        qd.r.wk <- qd.r.wk + 10^theta[i]*qd.r[,,i]
    }
    theta.old <- theta
    if (!is.null(random)) theta <- c(theta,zeta)
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
            warning("gss warning in sscox: CV iteration fails to converge")
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
    q.wk <- qd.r.wk <- 0
    for (i in 1:nq) {
        q.wk <- q.wk + 10^theta[i]*r[id.basis,,i]
        qd.r.wk <- qd.r.wk + 10^theta[i]*qd.r[,,i]
    }
    if (is.null(random)) q.wk0 <- 10^(lambda)*q.wk
    else {
        q.wk0 <- matrix(0,nxiz,nxiz)
        q.wk0[1:nxi,1:nxi] <- 10^(lambda)*q.wk
        q.wk0[(nxi+1):nxiz,(nxi+1):nxiz] <-
          10^(2*ran.scal)*random$sigma$fun(zeta,random$sigma$env)
        qd.r.wk <- cbind(qd.r.wk,10^ran.scal*random$qd.z)
    }
    se.aux <- .Fortran("coxaux",
                       as.double(cd), as.integer(nn),
                       as.double(q.wk0), as.integer(nxiz),
                       as.double(cbind(qd.r.wk,qd.s)), as.integer(nqd),
                       as.integer(bias$nt), as.double(bias$wt),
                       as.double(.Machine$double.eps),
                       as.double(qd.wt*bias$qd.wt),
                       double(nqd*bias$nt), double(bias$nt),
                       double(nn), v=double(nn*nn), double(nn*nn),
                       jpvt=integer(nn), PACKAGE="gss")[c("v","jpvt")]
    c <- cd[1:nxi]
    if (nz) b <- 10^ran.scal*cd[nxi+(1:nz)]
    else b <- NULL
    if (nnull) d <- cd[nxiz+(1:nnull)]
    else d <- NULL
    return(list(lambda=lambda,zeta=zeta,theta=theta,ran.scal=ran.scal,
                c=c,b=b,d=d,cv=cv,se.aux=se.aux))
}
