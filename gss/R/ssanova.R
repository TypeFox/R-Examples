## Fit ssanova model
ssanova <- function(formula,type=NULL,data=list(),weights,subset,
                    offset,na.action=na.omit,partial=NULL,
                    method="v",alpha=1.4,varht=1,
                    id.basis=NULL,nbasis=NULL,seed=NULL,random=NULL,
                    skip.iter=FALSE)
{
    ## Obtain model frame and model terms
    mf <- match.call()
    mf$type <- mf$method <- mf$varht <- mf$partial <- NULL
    mf$alpha <- mf$id.basis <- mf$nbasis <- mf$seed <- NULL
    mf$random <- mf$skip.iter <- NULL
    mf[[1]] <- as.name("model.frame")
    mf <- eval(mf,parent.frame())
    wt <- model.weights(mf)
    ## Generate sub-basis
    nobs <- dim(mf)[1]
    if (is.null(id.basis)) {
        if (is.null(nbasis))  nbasis <- max(30,ceiling(10*nobs^(2/9)))
        if (nbasis>=nobs)  nbasis <- nobs
        if (!is.null(seed))  set.seed(seed)
        id.basis <- sample(nobs,nbasis,prob=wt)
    }
    else {
        if (max(id.basis)>nobs|min(id.basis)<1)
            stop("gss error in ssanova: id.basis out of range")
        nbasis <- length(id.basis)
    }
    ## Generate terms
    term <- mkterm(mf,type)
    ## Generate random
    if (!is.null(random)) {
        if (class(random)=="formula") random <- mkran(random,data)
    }
    ## Generate s, r, and y
    s <- r <- NULL
    nq <- 0
    for (label in term$labels) {
        if (label=="1") {
            s <- cbind(s,rep(1,len=nobs))
            next
        }
        x <- mf[,term[[label]]$vlist]
        x.basis <- mf[id.basis,term[[label]]$vlist]
        nphi <- term[[label]]$nphi
        nrk <- term[[label]]$nrk
        if (nphi) {
            phi <- term[[label]]$phi
            for (i in 1:nphi)
                s <- cbind(s,phi$fun(x,nu=i,env=phi$env))
        }
        if (nrk) {
            rk <- term[[label]]$rk
            for (i in 1:nrk) {
                nq <- nq+1
                r <- array(c(r,rk$fun(x,x.basis,nu=i,env=rk$env,out=TRUE)),c(nobs,nbasis,nq))
            }
        }
    }
    if (is.null(r))
        stop("gss error in ssanova: use lm for models with only unpenalized terms")
    ## Add the partial term
    if (!is.null(partial)) {
        mf.p <- model.frame(partial,data)
        for (lab in colnames(mf.p)) mf[,lab] <- mf.p[,lab]
        mt.p <- attr(mf.p,"terms")
        lab.p <- labels(mt.p)
        matx.p <- model.matrix(mt.p,data)[,-1,drop=FALSE]
        if (dim(matx.p)[1]!=dim(mf)[1])
            stop("gss error in ssanova: partial data are of wrong size")
        matx.p <- scale(matx.p)
        center.p <- attr(matx.p,"scaled:center")
        scale.p <- attr(matx.p,"scaled:scale")
        s <- cbind(s,matx.p)
        part <- list(mt=mt.p,center=center.p,scale=scale.p)
    }
    else part <- lab.p <- NULL
    if (qr(s)$rank<dim(s)[2])
        stop("gss error in ssanova: unpenalized terms are linearly dependent")
    ## Prepare the data
    y <- model.response(mf,"numeric")
    offset <- model.offset(mf)
    if (!is.null(offset)) {
        term$labels <- c(term$labels,"offset")
        term$offset <- list(nphi=0,nrk=0)
        y <- y - offset
    }
    if (!is.null(wt)) wt <- sqrt(wt)
    ## Fit the model
    if (nq==1) {
        r <- r[,,1]
        z <- sspreg1(s,r,r[id.basis,],y,wt,method,alpha,varht,random)
    }
    else z <- mspreg1(s,r,id.basis,y,wt,method,alpha,varht,random,skip.iter)
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
    obj <- c(list(call=match.call(),mf=mf,terms=term,desc=desc,alpha=alpha,
                  id.basis=id.basis,partial=part,lab.p=lab.p,random=random,
                  skip.iter=skip.iter),z)
    class(obj) <- c("ssanova")
    obj
}

## Fit Single Smoothing Parameter (Gaussian) REGression
sspreg1 <- function(s,r,q,y,wt,method,alpha,varht,random)
{
    qr.trace <- FALSE
    if ((alpha<0)&(method%in%c("u","v"))) qr.trace <- TRUE
    alpha <- abs(alpha)
    ## get dimensions
    nobs <- nrow(r)
    nxi <- ncol(r)
    if (!is.null(s)) {
        if (is.vector(s)) nnull <- 1
        else nnull <- ncol(s)
    }
    else nnull <- 0
    if (!is.null(random)) nz <- ncol(as.matrix(random$z))
    else nz <- 0
    nxiz <- nxi + nz
    nn <- nxiz + nnull
    if (!is.null(wt)) {
        y <- wt*y
        s <- wt*s
        r <- wt*r
        if (!is.null(random)) random$z <- wt*random$z
    }
    ## cv function
    cv <- function(lambda) {
        if (is.null(random)) q.wk <- 10^(lambda+theta)*q
        else {
            q.wk <- matrix(0,nxiz,nxiz)
            q.wk[1:nxi,1:nxi] <- 10^(lambda[1]+theta)*q
            q.wk[(nxi+1):nxiz,(nxi+1):nxiz] <-
                10^(2*ran.scal)*random$sigma$fun(lambda[-1],random$sigma$env)
        }
        if (qr.trace) {
            qq.wk <- chol(q.wk,pivot=TRUE)
            sr <- cbind(s,10^theta*r[,attr(qq.wk,"pivot")])
            sr <- rbind(sr,cbind(matrix(0,nxiz,nnull),qq.wk))
            sr <- qr(sr,tol=0)
            rss <- mean(qr.resid(sr,c(y,rep(0,nxiz)))[1:nobs]^2)
            trc <- sum(qr.Q(sr)[1:nobs,]^2)/nobs
            if (method=="u") score <- rss + alpha*2*varht*trc
            if (method=="v") score <- rss/(1-alpha*trc)^2
            alpha.wk <- max(0,log.la0-lambda[1]-5)*(3-alpha) + alpha
            alpha.wk <- min(alpha.wk,3)
            if (alpha.wk>alpha) {
                if (method=="u") score <- score + (alpha.wk-alpha)*2*varht*trc
                if (method=="v") score <- rss/(1-alpha.wk*trc)^2
            }
            if (return.fit) {
                z <- .Fortran("reg",
                          as.double(cbind(s,10^theta*r)), as.integer(nobs), as.integer(nnull),
                          as.double(q.wk), as.integer(nxiz), as.double(y),
                          as.integer(switch(method,"u"=1,"v"=2,"m"=3)),
                          as.double(alpha), varht=as.double(varht),
                          score=double(1), dc=double(nn),
                          as.double(.Machine$double.eps),
                          chol=double(nn*nn), double(nn),
                          jpvt=as.integer(c(rep(1,nnull),rep(0,nxiz))),
                          wk=double(3*nobs+nnull+nz), rkv=integer(1), info=integer(1),
                          PACKAGE="gss")[c("score","varht","dc","chol","jpvt","wk","rkv","info")]
                z$score <- score
                assign("fit",z[c(1:5,7)],inherits=TRUE)
            }
        }
        else {
            z <- .Fortran("reg",
                          as.double(cbind(s,10^theta*r)), as.integer(nobs), as.integer(nnull),
                          as.double(q.wk), as.integer(nxiz), as.double(y),
                          as.integer(switch(method,"u"=1,"v"=2,"m"=3)),
                          as.double(alpha), varht=as.double(varht),
                          score=double(1), dc=double(nn),
                          as.double(.Machine$double.eps),
                          chol=double(nn*nn), double(nn),
                          jpvt=as.integer(c(rep(1,nnull),rep(0,nxiz))),
                          wk=double(3*nobs+nnull+nz), rkv=integer(1), info=integer(1),
                          PACKAGE="gss")[c("score","varht","dc","chol","jpvt","wk","rkv","info")]
            if (z$info) stop("gss error in ssanova: evaluation of GML score fails")
            assign("fit",z[c(1:5,7)],inherits=TRUE)
            score <- z$score
            alpha.wk <- max(0,log.la0-lambda[1]-5)*(3-alpha) + alpha
            alpha.wk <- min(alpha.wk,3)
            if (alpha.wk>alpha) {
                if (method=="u") score <- score + (alpha.wk-alpha)*2*varht*z$wk[2]
                if (method=="v") score <- z$wk[1]/(1-alpha.wk*z$wk[2])^2
            }
        }
        score
    }
    cv.wk <- function(lambda) cv.scale*cv(lambda)+cv.shift
    ## initialization
    tmp <- sum(r^2)
    if (is.null(s)) theta <- 0
    else theta <- log10(sum(s^2)/nnull/tmp*nxi) / 2
    log.la0 <- log10(tmp/sum(diag(q))) + theta
    if (!is.null(random)) {
        ran.scal <- theta - log10(sum(random$z^2)/nz/tmp*nxi) / 2
        r <- cbind(r,10^(ran.scal-theta)*random$z)
    }
    else ran.scal <- NULL
    ## lambda search
    return.fit <- FALSE
    fit <- NULL
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
                warning("gss warning in ssanova: iteration for model selection fails to converge")
                break
            }
        }
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
    }
    ## return
    return.fit <- TRUE
    jk1 <- cv(zz$est)
    if (is.null(random)) q.wk <- 10^theta*q
    else {
        q.wk <- matrix(0,nxiz,nxiz)
        q.wk[1:nxi,1:nxi] <- 10^theta*q
        q.wk[(nxi+1):nxiz,(nxi+1):nxiz] <-
            10^(2*ran.scal-zz$est[1])*random$sigma$fun(zz$est[-1],random$sigma$env)
    }
    se.aux <- regaux(s,10^theta*r,q.wk,zz$est[1],fit)
    c <- fit$dc[nnull+(1:nxi)]
    if (nnull) d <- fit$dc[1:nnull]
    else d <- NULL
    if (nz) b <- 10^(ran.scal)*fit$dc[nnull+nxi+(1:nz)]
    else b <- NULL
    c(list(method=method,theta=theta,ran.scal=ran.scal,c=c,d=d,b=b,
           nlambda=zz$est[1],zeta=zz$est[-1]),fit[-3],list(se.aux=se.aux))
}

## Fit Multiple Smoothing Parameter (Gaussian) REGression
mspreg1 <- function(s,r,id.basis,y,wt,method,alpha,varht,random,skip.iter)
{
    qr.trace <- FALSE
    if ((alpha<0)&(method%in%c("u","v"))) qr.trace <- TRUE
    alpha <- abs(alpha)
    ## get dimensions
    nobs <- nrow(r)
    nxi <- ncol(r)
    if (!is.null(s)) {
        if (is.vector(s)) nnull <- 1
        else nnull <- ncol(s)
    }
    else nnull <- 0
    if (!is.null(random)) nz <-ncol(as.matrix(random$z))
    else nz <- 0
    nxiz <- nxi + nz
    nn <- nxiz + nnull
    nq <- dim(r)[3]
    ## cv function
    cv <- function(theta) {
        ind.wk <- theta[1:nq]!=theta.old
        if (sum(ind.wk)==nq) {
            r.wk0 <- 0
            for (i in 1:nq) {
                r.wk0 <- r.wk0 + 10^theta[i]*r[,,i]
            }
            assign("r.wk",r.wk0+0,inherits=TRUE)
            assign("theta.old",theta[1:nq]+0,inherits=TRUE)
        }
        else {
            r.wk0 <- r.wk
            for (i in (1:nq)[ind.wk]) {
                theta.wk <- (10^(theta[i]-theta.old[i])-1)*10^theta.old[i]
                r.wk0 <- r.wk0 + theta.wk*r[,,i]
            }
        }
        qq.wk <- r.wk0[id.basis,]
        if (is.null(random)) q.wk <- 10^nlambda*qq.wk
        else {
            r.wk0 <- cbind(r.wk0,10^(ran.scal)*random$z)
            q.wk <- matrix(0,nxiz,nxiz)
            q.wk[1:nxi,1:nxi] <- 10^nlambda*qq.wk
            q.wk[(nxi+1):nxiz,(nxi+1):nxiz] <-
                10^(2*ran.scal)*random$sigma$fun(theta[-(1:nq)],random$sigma$env)
        }
        if (!is.null(wt)) {
            y.wk <- wt*y
            s.wk <- wt*s
            r.wk0 <- wt*r.wk0
        }
        if (qr.trace) {
            qq.wk <- chol(q.wk,pivot=TRUE)
            sr <- cbind(s.wk,r.wk0[,attr(qq.wk,"pivot")])
            sr <- rbind(sr,cbind(matrix(0,nxiz,nnull),qq.wk))
            sr <- qr(sr,tol=0)
            rss <- mean(qr.resid(sr,c(y.wk,rep(0,nxiz)))[1:nobs]^2)
            trc <- sum(qr.Q(sr)[1:nobs,]^2)/nobs
            if (method=="u") score <- rss + alpha*2*varht*trc
            if (method=="v") score <- rss/(1-alpha*trc)^2
            alpha.wk <- max(0,theta[1:nq]-log.th0-5)*(3-alpha) + alpha
            alpha.wk <- min(alpha.wk,3)
            if (alpha.wk>alpha) {
                if (method=="u") score <- score + (alpha.wk-alpha)*2*varht*trc
                if (method=="v") score <- rss/(1-alpha.wk*trc)^2
            }
            if (return.fit) {
                z <- .Fortran("reg",
                          as.double(cbind(s.wk,r.wk0)), as.integer(nobs), as.integer(nnull),
                          as.double(q.wk), as.integer(nxiz), as.double(y.wk),
                          as.integer(switch(method,"u"=1,"v"=2,"m"=3)),
                          as.double(alpha), varht=as.double(varht),
                          score=double(1), dc=double(nn),
                          as.double(.Machine$double.eps),
                          chol=double(nn*nn), double(nn),
                          jpvt=as.integer(c(rep(1,nnull),rep(0,nxiz))),
                          wk=double(3*nobs+nnull+nz), rkv=integer(1), info=integer(1),
                          PACKAGE="gss")[c("score","varht","dc","chol","jpvt","wk","rkv","info")]
                z$score <- score
                assign("fit",z[c(1:5,7)],inherits=TRUE)
            }
        }
        else {
            z <- .Fortran("reg",
                          as.double(cbind(s.wk,r.wk0)), as.integer(nobs), as.integer(nnull),
                          as.double(q.wk), as.integer(nxiz), as.double(y.wk),
                          as.integer(switch(method,"u"=1,"v"=2,"m"=3)),
                          as.double(alpha), varht=as.double(varht),
                          score=double(1), dc=double(nn),
                          as.double(.Machine$double.eps),
                          chol=double(nn*nn), double(nn),
                          jpvt=as.integer(c(rep(1,nnull),rep(0,nxiz))),
                          wk=double(3*nobs+nnull+nz), rkv=integer(1), info=integer(1),
                          PACKAGE="gss")[c("score","varht","dc","chol","jpvt","wk","rkv","info")]
            if (z$info) stop("gss error in ssanova: evaluation of GML score fails")
            assign("fit",z[c(1:5,7)],inherits=TRUE)
            score <- z$score
            alpha.wk <- max(0,theta[1:nq]-log.th0-5)*(3-alpha) + alpha
            alpha.wk <- min(alpha.wk,3)
            if (alpha.wk>alpha) {
                if (method=="u") score <- score + (alpha.wk-alpha)*2*varht*z$wk[2]
                if (method=="v") score <- z$wk[1]/(1-alpha.wk*z$wk[2])^2
            }
        }
        score
    }
    cv.wk <- function(theta) cv.scale*cv(theta)+cv.shift
    ## initialization
    theta <- -log10(apply(r[id.basis,,],3,function(x)sum(diag(x))))
    r.wk <- 0
    for (i in 1:nq) {
        r.wk <- r.wk + 10^theta[i]*r[,,i]
    }
    ## theta adjustment
    return.fit <- FALSE
    z <- sspreg1(s,r.wk,r.wk[id.basis,],y,wt,method,alpha,varht,random)
    theta <- theta + z$theta
    r.wk <- 0
    for (i in 1:nq) {
        theta[i] <- 2*theta[i] + log10(t(z$c)%*%r[id.basis,,i]%*%z$c)
        r.wk <- r.wk + 10^theta[i]*r[,,i]
    }
    if (!is.null(wt)) q.wk <- wt*r.wk
    else q.wk <- r.wk
    log.la0 <- log10(sum(q.wk^2)/sum(diag(r.wk[id.basis,])))
    log.th0 <- theta-log.la0
    ## lambda search
    z <- sspreg1(s,r.wk,r.wk[id.basis,],y,wt,method,alpha,varht,random)
    nlambda <- z$nlambda
    log.th0 <- log.th0 + z$nlambda
    theta <- theta + z$theta
    if (!is.null(random)) ran.scal <- z$ran.scal
    ## early return
    if (skip.iter) {
        z$theta <- theta
        return(z)
    }
    ## theta search
    fit <- NULL
    counter <- 0
    y.wk <- y
    s.wk <- s
    r.wk <- 0
    for (i in 1:nq) {
        r.wk <- r.wk + 10^theta[i]*r[,,i]
    }
    theta.old <- theta
    if (!is.null(random)) theta <- c(theta,z$zeta)
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
            warning("gss warning in ssanova: iteration for model selection fails to converge")
            break
        }
    }
    ## return
    return.fit <- TRUE
    jk1 <- cv(zz$est)
    r.wk <- 0
    for (i in 1:nq) {
        r.wk <- r.wk + 10^zz$est[i]*r[,,i]
    }
    qq.wk <- r.wk[id.basis,]
    if (is.null(random)) q.wk <- qq.wk
    else {
        r.wk <- cbind(r.wk,10^(ran.scal)*random$z)
        q.wk <- matrix(0,nxiz,nxiz)
        q.wk[1:nxi,1:nxi] <- qq.wk
        q.wk[(nxi+1):nxiz,(nxi+1):nxiz] <-
            10^(2*ran.scal-nlambda)*random$sigma$fun(zz$est[-(1:nq)],random$sigma$env)
    }
    if (!is.null(wt)) {
        s <- wt*s
        r.wk <- wt*r.wk
    }
    se.aux <- regaux(s,r.wk,q.wk,nlambda,fit)
    c <- fit$dc[nnull+(1:nxi)]
    if (nnull) d <- fit$dc[1:nnull]
    else d <- NULL
    if (nz) b <- 10^(ran.scal)*fit$dc[nnull+nxi+(1:nz)]
    else b <- NULL
    c(list(method=method,theta=zz$est[1:nq],c=c,d=d,b=b,nlambda=nlambda,
           zeta=zz$est[-(1:nq)]),fit[-3],list(se.aux=se.aux))
}

## Auxiliary Quantities for Standard Error Calculation
regaux <- function(s,r,q,nlambda,fit)
{
    nnull <- dim(s)[2]
    nn <- nnull +  dim(q)[1]
    zzz <- eigen(q,symmetric=TRUE)
    rkq <- min(fit$rkv-nnull,sum(zzz$val/zzz$val[1]>sqrt(.Machine$double.eps)))
    val <- zzz$val[1:rkq]
    vec <- zzz$vec[,1:rkq,drop=FALSE]
    if (nnull) {
        wk1 <- qr(s)
        wk1 <- (qr.qty(wk1,r%*%vec))[-(1:nnull),]
    }
    else wk1 <- r%*%vec
    wk2 <- t(t(wk1)/sqrt(val))
    wk2 <- t(wk2)%*%wk2
    wk2 <- solve(wk2+diag(10^nlambda,dim(wk2)[1]),wk2)
    wk2 <- (wk2+t(wk2))/2
    wk2 <- t(wk2/sqrt(val))/sqrt(val)
    wk2 <- diag(1/val,dim(wk2)[1])-wk2
    z <- .Fortran("regaux",
                  as.double(fit$chol), as.integer(nn),
                  as.integer(fit$jpvt), as.integer(fit$rkv),
                  drcr=as.double(t(cbind(s,r))%*%r%*%vec), as.integer(rkq),
                  sms=double(nnull^2), as.integer(nnull), double(nn*nnull),
                  PACKAGE="gss")[c("drcr","sms")]
    drcr <- matrix(z$drcr,nn,rkq)
    dr <- drcr[1:nnull,,drop=FALSE]
    sms <- 10^nlambda*matrix(z$sms,nnull,nnull)
    wk1 <- matrix(0,nnull+rkq,nnull+rkq)
    wk1[1:nnull,1:nnull] <- sms
    wk1[1:nnull,nnull+(1:rkq)] <- -t(t(dr)/val)
    wk1[nnull+(1:rkq),nnull+(1:rkq)] <- wk2
    z <- chol(wk1,pivot=TRUE)
    wk1 <- z
    rkw <- attr(z,"rank")
    while (wk1[rkw,rkw]<wk1[1,1]*sqrt(.Machine$double.eps)) rkw <- rkw-1
    wk1[row(wk1)>col(wk1)] <- 0
    if (rkw<nnull+rkq)
        wk1[(rkw+1):(nnull+rkq),(rkw+1):(nnull+rkq)] <- diag(0,nnull+rkq-rkw)
    hfac <- wk1
    hfac[,attr(z,"pivot")] <- wk1
    list(vec=vec,hfac=hfac)
}
