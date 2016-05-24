## Fit hazard model
sshzd1 <- function(formula,type=NULL,data=list(),alpha=1.4,
                   weights=NULL,subset,na.action=na.omit,rho="marginal",
                   partial=NULL,id.basis=NULL,nbasis=NULL,seed=NULL,
                   random=NULL,prec=1e-7,maxiter=30,skip.iter=FALSE)
{
    ## Local functions handling formula
    Surv <- function(time,status,start=0) {
        tname <- as.character(as.list(match.call())$time)
        if (!is.numeric(time)|!is.vector(time))
            stop("gss error in sshzd1: time should be a numerical vector")
        if ((nobs <- length(time))-length(status))
            stop("gss error in sshzd1: time and status mismatch in size")
        if ((length(start)-nobs)&(length(start)-1))
            stop("gss error in sshzd1: time and start mismatch in size")
        if (any(start>time))
            stop("gss error in sshzd1: start after follow-up time")
        if (min(start)<0)
            stop("gss error in sshzd1: start before time 0")
        time <- cbind(start,time)
        list(tname=tname,start=time[,1],end=time[,2],status=as.logical(status))
    }
    ## Obtain model frame and model terms
    mf <- match.call()
    mf$type <- mf$alpha <- mf$random <- mf$partial <- NULL
    mf$id.basis <- mf$nbasis <- mf$seed <- mf$rho <- NULL
    mf$prec <- mf$maxiter <- mf$skip.iter <- NULL
    term.wk <- terms.formula(formula)
    ## response
    resp <- attr(term.wk,"variable")[[2]]
    ind.wk <- length(strsplit(deparse(resp),'')[[1]])
    if ((substr(deparse(resp),1,5)!='Surv(')
        |(substr(deparse(resp),ind.wk,ind.wk)!=')'))
        stop("gss error in sshzd1: response should be Surv(...)")
    yy <- with(data,eval(resp))
    tname <- yy$tname
    ## model frame
    term.labels <- attr(term.wk,"term.labels")
    if (!(tname%in%term.labels))
        stop("gss error in sshzd1: time main effect missing in model")
    mf[[1]] <- as.name("model.frame")
    mf[[2]] <- eval(parse(text=paste("~",paste(term.labels,collapse="+"))))
    mf <- eval(mf,parent.frame())
    ## Use sshzd in lack of covariate
    if (all(tname==names(mf))) stop("use sshzd when covariate is absent")
    ## trim yy if subset is used
    nobs <- nrow(mf)
    if (nobs<length(yy$end)) {
        yy$start <- yy$start[subset]
        yy$end <- yy$end[subset]
        yy$status <- yy$status[subset]
    }
    ## Generate sub-basis
    cnt <- model.weights(mf)
    if (is.null(cnt)) yy$cnt <- rep(1,nobs)
    else {
        yy$cnt <- cnt
        mf["(weights)"] <- NULL
    }
    if (is.null(id.basis)) {
        if (is.null(nbasis)) nbasis <- max(30,ceiling(10*nobs^(2/9)))
        if (nbasis>sum(yy$status)) nbasis <- sum(yy$status)
        if (!is.null(seed)) set.seed(seed)
        id.basis <- sample((1:nobs)[yy$status],nbasis,prob=cnt[yy$status])
    }
    else {
        if (!all(id.basis%in%(1:nobs)[yy$status]))
            stop("gss error in sshzd1: id.basis not all at failure cases")
        nbasis <- length(id.basis)
    }
    id.wk <- NULL
    nT <- sum(yy$status)
    for (i in 1:nbasis) {
        id.wk <- c(id.wk,(1:nT)[(1:nobs)[yy$status]%in%id.basis[i]])
    }
    ## set domain and type for time
    mn <- min(yy$start)
    mx <- max(yy$end)
    tdomain <- c(max(mn-.05*(mx-mn),0),mx)
    if (is.null(type[[tname]])) type[[tname]] <- list("cubic",tdomain)
    if (length(type[[tname]])==1) type[[tname]] <- c(type[[tname]],tdomain)
    if (!(type[[tname]][[1]]%in%c("cubic","linear")))
        stop("gss error in sshzd1: wrong type")
    if ((min(type[[tname]][[2]])>min(tdomain))|
        (max(type[[tname]][[2]])<max(tdomain)))
        stop("gss error in sshzd1: time range not covering domain")
    ## Generate terms    
    term <- mkterm(mf,type)
    ## Generate random
    if (!is.null(random)) {
        if (class(random)=="formula") random <- mkran(random,data)
    }
    ## Generate Gauss-Legendre quadrature
    nmesh <- 200
    quad <- gauss.quad(nmesh,tdomain)
    ## Generate partial terms
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
        part <- list(mt=mt.p,center=center.p,scale=scale.p)
    }
    else part <- lab.p <- NULL
    ## obtain unique covariate observations
    xnames <- names(mf)
    xnames <- xnames[!xnames%in%tname]
    if (length(xnames)) {
        xx <- mf[,xnames,drop=FALSE]
        if (!is.null(partial)) xx <- cbind(xx,matx.p)
        if (!is.null(random)) xx <- cbind(xx,random$z)
        xx <- apply(xx,1,function(x)paste(x,collapse="\r"))
        ux <- unique(xx)
        nx <- length(ux)
        x.dup.ind <- duplicated(xx)
        x.dup <- as.vector(xx[x.dup.ind])
        x.pt <- mf[!x.dup.ind,xnames,drop=FALSE]
        ## xx[i,]==x.pt[x.ind[i],]
        x.ind <- 1:nobs
        x.ind[!x.dup.ind] <- 1:nx
        if (nobs-nx) {
            x.ind.wk <- range <- 1:(nobs-nx)
            for (i in 1:nx) {
                range.wk <- NULL
                for (j in range) {
                    if (identical(ux[i],x.dup[j])) {
                        x.ind.wk[j] <- i
                        range.wk <- c(range.wk,j)
                    }
                }
                if (!is.null(range.wk)) range <- range[!(range%in%range.wk)]
            }
            x.ind[x.dup.ind] <- x.ind.wk
        }
        if (!is.null(random)) {
            qd.z <- random$z[!x.dup.ind,]
            random$z <- random$z[yy$status,]
        }
    }
    else stop("gss error in sshzd1: missing covariate")
    ## calculate rho
    if (rho=="marginal") {
        rho.wk <- sshzd(Surv(end,status,start)~end,data=yy,
                        id.basis=id.basis,weights=cnt,alpha=2)
        rho.qd <- hzdcurve.sshzd(rho.wk,quad$pt)
        rhowk <- hzdcurve.sshzd(rho.wk,yy$end[yy$status])
    }
    if (rho=="weibull") {
        y.wk <- cbind(yy$end,yy$status,yy$start)
        form <- as.formula(paste("y.wk~",paste(xnames,collapse="+")))
        rho.wk <- gssanova(form,family="weibull",partial=partial,data=data,
                           id.basis=id.basis,weights=cnt,alpha=2)
        yhat <- predict(rho.wk,rho.wk$mf)
        rho.qd <- exp(rho.wk$nu*outer(log(quad$pt),yhat[!x.dup.ind],"-"))/quad$pt
        rhowk <- (exp(rho.wk$nu*(log(yy$end)-yhat))/yy$end)[yy$status]
    }
    ## integration weights at x.pt[i,]
    qd.wt <- matrix(0,nmesh,nx)
    for (i in 1:nobs) {
        wk <- (quad$pt<=yy$end[i])&(quad$pt>yy$start[i])
        if (is.vector(rho.qd)) wk <- wk*rho.qd
        else wk <- wk*rho.qd[,x.ind[i]]
        if (is.null(cnt)) qd.wt[,x.ind[i]] <- qd.wt[,x.ind[i]]+wk
        else qd.wt[,x.ind[i]] <- qd.wt[,x.ind[i]]+cnt[i]*wk
    }
    if (is.null(cnt)) qd.wt <- quad$wt*qd.wt/nobs
    else qd.wt <- quad$wt*qd.wt/sum(cnt)
    ## Generate s, r, int.s, and int.r
    s <- r <- int.s <- int.r <- NULL
    nq <- 0
    for (label in term$labels) {
        if (label=="1") {
            s <- cbind(s,rep(1,len=nT))
            int.s <- c(int.s,sum(qd.wt))
            next
        }
        vlist <- term[[label]]$vlist
        x.list <- xnames[xnames%in%vlist]
        xy <- mf[yy$status,vlist]
        xy.basis <- mf[id.basis,vlist]
        qd.xy <- data.frame(matrix(0,nmesh,length(vlist)))
        names(qd.xy) <- vlist
        if (tname%in%vlist) qd.xy[,tname] <- quad$pt
        if (length(x.list)) xx <- x.pt[,x.list,drop=FALSE]
        else xx <- NULL
        nphi <- term[[label]]$nphi
        nrk <- term[[label]]$nrk
        if (nphi) {
            phi <- term[[label]]$phi
            for (i in 1:nphi) {
                s <- cbind(s,phi$fun(xy,nu=i,env=phi$env))
                if (is.null(xx)) {
                    qd.wk <- phi$fun(qd.xy[,,drop=TRUE],nu=i,env=phi$env)
                    int.s <- c(int.s,sum(qd.wk*apply(qd.wt,1,sum)))
                }
                else {
                    int.s.wk <- 0
                    for (j in 1:nx) {
                        qd.xy[,x.list] <- xx[rep(j,nmesh),]
                        qd.wk <- phi$fun(qd.xy[,,drop=TRUE],i,phi$env)
                        int.s.wk <- int.s.wk + sum(qd.wk*qd.wt[,j])
                    }
                    int.s <- c(int.s,int.s.wk)
                }
            }
        }
        if (nrk) {
            rk <- term[[label]]$rk
            for (i in 1:nrk) {
                nq <- nq+1
                r <- array(c(r,rk$fun(xy,xy.basis,nu=i,env=rk$env,out=TRUE)),c(nT,nbasis,nq))
                if (is.null(xx)) {
                    qd.wk <- rk$fun(qd.xy[,,drop=TRUE],xy.basis,i,rk$env,out=TRUE)
                    int.r <- cbind(int.r,apply(apply(qd.wt,1,sum)*qd.wk,2,sum))
                }
                else {
                    int.r.wk <- 0
                    for (j in 1:nx) {
                        qd.xy[,x.list] <- xx[rep(j,nmesh),]
                        qd.wk <- rk$fun(qd.xy[,,drop=TRUE],xy.basis,i,rk$env,TRUE)
                        int.r.wk <- int.r.wk + apply(qd.wt[,j]*qd.wk,2,sum)
                    }
                    int.r <- cbind(int.r,int.r.wk)
                }
            }
        }
    }
    ## Add the partial term
    if (!is.null(partial)) {
        s <- cbind(s,matx.p[yy$status,])
        int.s <- c(int.s,t(matx.p[!x.dup.ind,])%*%apply(qd.wt,2,sum))
        part$pt <- matx.p[!x.dup.ind,,drop=FALSE]
    }
    ## generate int.z
    if (!is.null(random)) random$int.z <- t(qd.z)%*%apply(qd.wt,2,sum)
    ## Check s rank
    if (!is.null(s)) {
        nnull <- dim(s)[2]
        if (qr(s)$rank<nnull)
            stop("gss error in sshzd1: unpenalized terms are linearly dependent")
    }
    ## Fit the model
    Nobs <- ifelse(is.null(cnt),nobs,sum(cnt))
    if (!is.null(cnt)) cntt <- cnt[yy$status]
    else cntt <- NULL
    z <- msphzd1(s,r,id.wk,Nobs,cntt,int.s,int.r,rhowk,random,prec,maxiter,alpha,skip.iter)
    ## cfit
    if (!is.null(random)) rhowk <- rhowk*exp(-random$z%*%z$b)
    if (!is.null(cnt)) cfit <- sum(cntt*rhowk)/Nobs/sum(qd.wt)
    else cfit <- sum(rhowk)/Nobs/sum(qd.wt)
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
    obj <- c(list(call=match.call(),mf=mf,cnt=cnt,terms=term,desc=desc,yy=yy,
                  alpha=alpha,tname=tname,xnames=xnames,tdomain=tdomain,cfit=cfit,
                  quad=quad,x.pt=x.pt,qd.wt=qd.wt,id.basis=id.basis,partial=part,
                  lab.p=lab.p,random=random,int.s=int.s,int.r=int.r,skip.iter=skip.iter),z)
    obj$se.aux$v <- sqrt(Nobs)*obj$se.aux$v
    class(obj) <- c("sshzd1","sshzd")
    obj
}

## Fit (multiple smoothing parameter) hazard function
msphzd1 <- function(s,r,id.wk,Nobs,cnt,int.s,int.r,rho,random,prec,maxiter,alpha,skip.iter)
{
    nT <- dim(r)[1]
    nxi <- dim(r)[2]
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
        fit <- .Fortran("hzdnewton10",
                        cd=as.double(cd), as.integer(nn),
                        as.double(q.wk0), as.integer(nxiz),
                        as.double(cbind(r.wk,s)), as.integer(nT),
                        as.integer(Nobs), as.integer(sum(cnt)),
                        as.integer(cnt),
                        as.double(c(int.r.wk,int.s)), as.double(rho),
                        as.double(prec), as.integer(maxiter),
                        as.double(.Machine$double.eps), integer(nn),
                        wk=double(2*nT+nn*(nn+3)),
                        info=integer(1),PACKAGE="gss")
        if (fit$info==1) stop("gss error in sshzd1: Newton iteration diverges")
        if (fit$info==2) warning("gss warning in sshzd1: Newton iteration fails to converge")
        assign("cd",fit$cd,inherits=TRUE)
        cv <- fit$wk[1]+alpha*fit$wk[2]
        alpha.wk <- max(0,log.la0-lambda[1]-5)*(3-alpha) + alpha
        alpha.wk <- min(alpha.wk,3)
        adj <- ifelse (alpha.wk>alpha,(alpha.wk-alpha)*fit$wk[2],0)
        cv+adj
    }
    cv.s.wk <- function(lambda) cv.scale*cv.s(lambda)+cv.shift
    cv.m <- function(theta) {
        ind.wk <- theta[1:nq]!=theta.old
        if (sum(ind.wk)==nq) {
            r.wk0 <- int.r.wk0 <- 0
            for (i in 1:nq) {
                r.wk0 <- r.wk0 + 10^theta[i]*r[,,i]
                int.r.wk0 <- int.r.wk0 + 10^theta[i]*int.r[,i]
            }
            assign("r.wk",r.wk0+0,inherits=TRUE)
            assign("int.r.wk",int.r.wk0+0,inherits=TRUE)
            assign("theta.old",theta[1:nq]+0,inherits=TRUE)
        }
        else {
            r.wk0 <- r.wk
            int.r.wk0 <- int.r.wk
            for (i in (1:nq)[ind.wk]) {
                theta.wk <- (10^(theta[i]-theta.old[i])-1)*10^theta.old[i]
                r.wk0 <- r.wk0 + theta.wk*r[,,i]
                int.r.wk0 <- int.r.wk0 + theta.wk*int.r[,i]
            }
        }
        q.wk <- r.wk0[id.wk,]
        if (is.null(random)) q.wk0 <- 10^(lambda)*q.wk
        else {
            r.wk0 <- cbind(r.wk0,10^ran.scal*random$z)
            int.r.wk0 <- c(int.r.wk0,10^ran.scal*random$int.z)
            q.wk0 <- matrix(0,nxiz,nxiz)
            q.wk0[1:nxi,1:nxi] <- 10^(lambda)*q.wk
            q.wk0[(nxi+1):nxiz,(nxi+1):nxiz] <-
              10^(2*ran.scal)*random$sigma$fun(theta[-(1:nq)],random$sigma$env)
        }
        fit <- .Fortran("hzdnewton10",
                        cd=as.double(cd), as.integer(nn),
                        as.double(q.wk0), as.integer(nxiz),
                        as.double(cbind(r.wk0,s)), as.integer(nT),
                        as.integer(Nobs), as.integer(sum(cnt)),
                        as.integer(cnt),
                        as.double(c(int.r.wk0,int.s)), as.double(rho),
                        as.double(prec), as.integer(maxiter),
                        as.double(.Machine$double.eps), integer(nn),
                        wk=double(2*nT+nn*(nn+3)),
                        info=integer(1),PACKAGE="gss")
        if (fit$info==1) stop("gss error in sshzd: Newton iteration diverges")
        if (fit$info==2) warning("gss warning in sshzd: Newton iteration fails to converge")
        assign("cd",fit$cd,inherits=TRUE)
        cv <- fit$wk[1]+alpha*fit$wk[2]
        alpha.wk <- max(0,theta[1:nq]-log.th0-5)*(3-alpha) + alpha
        alpha.wk <- min(alpha.wk,3)
        adj <- ifelse (alpha.wk>alpha,(alpha.wk-alpha)*fit$wk[2],0)
        cv+adj
    }
    cv.m.wk <- function(theta) cv.scale*cv.m(theta)+cv.shift
    ## Initialization
    theta <- -log10(apply(r[id.wk,,,drop=FALSE],3,function(x)sum(diag(x))))
    nq <- length(theta)
    r.wk <- int.r.wk <- 0
    for (i in 1:nq) {
        r.wk <- r.wk + 10^theta[i]*r[,,i]
        int.r.wk <- int.r.wk + 10^theta[i]*int.r[,i]
    }
    v.r <- sum(rho*r.wk^2)
    if (nnull) {
        v.s <- sum(rho*s^2)
        theta.wk <- log10(v.s/nnull/v.r*nxi) / 2
    }
    else theta.wk <- 0
    if (!is.null(random)) {
        v.z <- sum(rho*random$z^2)
        ran.scal <- theta.wk - log10(v.z/nz/v.r*nxi) / 2
    }
    else ran.scal <- NULL
    theta <- theta + theta.wk
    r.wk <- 10^theta.wk*r.wk
    int.r.wk <- 10^theta.wk*int.r.wk
    q.wk <- r.wk[id.wk,]
    if (!is.null(random)) {
        r.wk <- cbind(r.wk,10^ran.scal*random$z)
        int.r.wk <- c(int.r.wk,10^ran.scal*random$int.z)
    }
    log.la0 <- log10(sum(v.r)/sum(diag(q.wk))) + 2*theta.wk
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
                warning("gss warning in sshzd: iteration for model selection fails to converge")
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
        if (!is.null(cnt)) rho <- rho*cnt
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
        se.aux <- .Fortran("hzdaux101",
                           as.double(cd), as.integer(nn),
                           as.double(q.wk0), as.integer(nxiz),
                           as.double(cbind(r.wk,s)), as.integer(nT),
                           as.double(rho/Nobs), as.double(.Machine$double.eps),
                           v=double(nn*nn), jpvt=integer(nn),
                           PACKAGE="gss")[c("v","jpvt")]
        c <- cd[1:nxi]
        if (nz) b <- 10^ran.scal*cd[nxi+(1:nz)]
        else b <- NULL
        if (nnull) d <- cd[nxiz+(1:nnull)]
        else d <- NULL
        return(list(lambda=lambda,zeta=zeta,theta=theta,ran.scal=ran.scal,
                    c=c,b=b,d=d,cv=cv,se.aux=se.aux))
    }
    ## theta adjustment
    for (i in 1:nq) {
        theta[i] <- 2*theta[i] + log10(t(cd[1:nxi])%*%r[id.wk,,i]%*%cd[1:nxi])
    }
    r.wk <- int.r.wk <- 0
    for (i in 1:nq) {
        r.wk <- r.wk + 10^theta[i]*r[,,i]
        int.r.wk <- int.r.wk + 10^theta[i]*int.r[,i]
    }
    v.r <- sum(rho*r.wk^2)
    if (nnull) {
        v.s <- sum(rho*s^2)
        theta.wk <- log10(v.s/nnull/v.r*nxi) / 2
    }
    else theta.wk <- 0
    if (!is.null(random)) {
        v.z <- sum(rho*random$z^2)
        ran.scal <- theta.wk - log10(v.z/nz/v.r*nxi) / 2
    }
    else ran.scal <- NULL
    theta <- theta + theta.wk
    r.wk <- 10^theta.wk*r.wk
    int.r.wk <- 10^theta.wk*int.r.wk
    q.wk <- r.wk[id.wk,]
    if (!is.null(random)) {
        r.wk <- cbind(r.wk,10^ran.scal*random$z)
        int.r.wk <- c(int.r.wk,10^ran.scal*random$int.z)
    }
    log.la0 <- log10(sum(v.r)/sum(diag(q.wk))) + 2*theta.wk
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
                warning("gss warning in sshzd: iteration for model selection fails to converge")
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
        if (!is.null(cnt)) rho <- rho*cnt
        if (is.null(random)) q.wk0 <- 10^(lambda)*q.wk
        else {
            q.wk0 <- matrix(0,nxiz,nxiz)
            q.wk0[1:nxi,1:nxi] <- 10^(lambda)*q.wk
            q.wk0[(nxi+1):nxiz,(nxi+1):nxiz] <-
                10^(2*ran.scal)*random$sigma$fun(zeta,random$sigma$env)
        }
        se.aux <- .Fortran("hzdaux101",
                           as.double(cd), as.integer(nn),
                           as.double(q.wk0), as.integer(nxiz),
                           as.double(cbind(r.wk,s)), as.integer(nT),
                           as.double(rho/Nobs), as.double(.Machine$double.eps),
                           v=double(nn*nn), jpvt=integer(nn),
                           PACKAGE="gss")[c("v","jpvt")]
        c <- cd[1:nxi]
        if (nz) b <- 10^ran.scal*cd[nxi+(1:nz)]
        else b <- NULL
        if (nnull) d <- cd[nxiz+(1:nnull)]
        else d <- NULL
        return(list(lambda=lambda,zeta=zeta,theta=theta,ran.scal=ran.scal,
                    c=c,b=b,d=d,cv=cv,se.aux=se.aux))
    }
    ## theta search
    counter <- 0
    r.wk <- int.r.wk <- 0
    for (i in 1:nq) {
        r.wk <- r.wk + 10^theta[i]*r[,,i]
        int.r.wk <- int.r.wk + 10^theta[i]*int.r[,i]
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
            warning("gss warning in sshzd1: CV iteration fails to converge")
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
    if (!is.null(cnt)) rho <- rho*cnt
    r.wk <- 0
    for (i in 1:nq) {
        r.wk <- r.wk + 10^theta[i]*r[,,i]
    }
    q.wk <- r.wk[id.wk,]
    if (is.null(random)) q.wk0 <- 10^(lambda)*q.wk
    else {
        r.wk <- cbind(r.wk,10^ran.scal*random$z)
        q.wk0 <- matrix(0,nxiz,nxiz)
        q.wk0[1:nxi,1:nxi] <- 10^(lambda)*q.wk
        q.wk0[(nxi+1):nxiz,(nxi+1):nxiz] <-
          10^(2*ran.scal)*random$sigma$fun(zeta,random$sigma$env)
    }
    se.aux <- .Fortran("hzdaux101",
                       as.double(cd), as.integer(nn),
                       as.double(q.wk0), as.integer(nxiz),
                       as.double(cbind(r.wk,s)), as.integer(nT),
                       as.double(rho/Nobs), as.double(.Machine$double.eps),
                       v=double(nn*nn), jpvt=integer(nn),
                       PACKAGE="gss")[c("v","jpvt")]
    c <- cd[1:nxi]
    if (nz) b <- 10^ran.scal*cd[nxi+(1:nz)]
    else b <- NULL
    if (nnull) d <- cd[nxiz+(1:nnull)]
    else d <- NULL
    list(lambda=lambda,zeta=zeta,theta=theta,ran.scal=ran.scal,
         c=c,b=b,d=d,cv=cv,se.aux=se.aux)
}
