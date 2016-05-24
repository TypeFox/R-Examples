## Fit hazard model
sshzd2d1 <- function(formula1,formula2,symmetry=FALSE,data,alpha=1.4,weights=NULL,
                     subset=NULL,rho="marginal",id.basis=NULL,nbasis=NULL,
                     seed=NULL,prec=1e-7,maxiter=30,skip.iter=FALSE)
{
    ## Prepare data
    if (!is.null(subset)) {
        data <- data[subset,]
        if (!is.null(weights)) {
            id.wk <- apply(!is.na(data),1,all)
            cnt <- weights[id.wk]
        }
        else cnt <- NULL
        data <- na.omit(data[subset,])
    }
    else {
        if (!is.null(weights)) {
            id.wk <- apply(!is.na(data),1,all)
            cnt <- weights[id.wk]
        }
        else cnt <- NULL
        data <- na.omit(data)
    }
    ## Extract formulas
    if (class(formula1)=="formula") {
        form1 <- formula1
        part1 <- random1 <- type1 <- NULL
    }
    else {
        if (class(formula1)!="list")
            stop("gss error in sshzd2d1: models must be specified via formulas or lists")
        form1 <- formula1[[1]]
        part1 <- formula1$partial
        random1 <- formula1$random
        type1 <- formula1$type
    }
    if (class(formula2)=="formula") {
        form2 <- formula2
        part2 <- random2 <- type2 <- NULL
    }
    else {
        if (class(formula2)!="list")
            stop("gss error in sshzd2d1: models must be specified via formulas or lists")
        form2 <- formula2[[1]]
        part2 <- formula2$partial
        random2 <- formula2$random
        type2 <- formula2$type
    }
    ## Local function handling formula
    Surv <- function(time,status,start=0) {
        tname <- as.character(as.list(match.call())$time)
        if (!is.numeric(time)|!is.vector(time))
            stop("gss error in sshzd2d1: time should be a numerical vector")
        if ((nobs <- length(time))-length(status))
            stop("gss error in sshzd2d1: time and status mismatch in size")
        if ((length(start)-nobs)&(length(start)-1))
            stop("gss error in sshzd2d1: time and start mismatch in size")
        if (any(start>time))
            stop("gss error in sshzd2d1: start after follow-up time")
        if (min(start)<0)
            stop("gss error in sshzd2d1: start before time 0")
        time <- cbind(start,time)
        list(tname=tname,start=time[,1],end=time[,2],status=as.logical(status))
    }
    ## Obtain model terms
    for (i in 1:2) {
        ## Formula
        form.wk <- switch(i,form1,form2)
        term.wk <- terms.formula(form.wk)
        resp <- attr(term.wk,"variable")[[2]]
        ind.wk <- length(strsplit(deparse(resp),'')[[1]])
        if ((substr(deparse(resp),1,5)!='Surv(')
            |(substr(deparse(resp),ind.wk,ind.wk)!=')'))
            stop("gss error in sshzd2d1: response should be Surv(...)")
        yy <- with(data,eval(resp))
        tname <- yy$tname
        term.labels <- attr(term.wk,"term.labels")
        if (!(tname%in%term.labels))
            stop("gss error in sshzd2d1: time main effect missing in model")
        form.wk <- eval(parse(text=paste("~",paste(term.labels,collapse="+"))))
        mf.wk <- model.frame(form.wk,data)
        ## Partial
        part.wk <- switch(i,part1,part2)
        if (!is.null(part.wk)) {
            mf.p.wk <- model.frame(part.wk,data)
            mt.p.wk <- attr(mf.p.wk,"terms")
            matx.p.wk <- model.matrix(mt.p.wk,data)[,-1,drop=FALSE]
            if (dim(matx.p.wk)[1]!=dim(mf.wk)[1])
                stop("gss error in sshzd2d1: partial data are of wrong size")
        }
        else mf.p.wk <- mt.p.wk <- matx.p.wk <- NULL
        ## Random
        random.wk <- switch(i,random1,random2)
        if (!is.null(random.wk)) {
            if (class(random.wk)=="formula") random.wk <- mkran(random.wk,data)
        }
        else random.wk <- NULL
        ## Set domain and type for time
        type.wk <- switch(i,type1,type2)
        mn <- min(yy$start)
        mx <- max(yy$end)
        tdomain <- c(max(mn-.05*(mx-mn),0),mx)
        tname <- yy$tname
        if (is.null(type.wk[[tname]])) type.wk[[tname]] <- list("cubic",tdomain)
        if (length(type.wk[[tname]])==1) type.wk[[tname]] <- c(type.wk[[tname]],tdomain)
        if (!(type.wk[[tname]][[1]]%in%c("cubic","linear")))
            stop("gss error in sshzd2d1: wrong type")
        if ((min(type.wk[[tname]][[2]])>min(tdomain))|
            (max(type.wk[[tname]][[2]])<max(tdomain)))
            stop("gss error in sshzd2d1: time range not covering domain")
        ## Save
        if (i==1) {
            mf1 <- mf.wk
            yy1 <- yy
            mf.p1 <- mf.p.wk
            mt.p1 <- mt.p.wk
            matx.p1 <- matx.p.wk
            random1 <- random.wk
            type1 <- type.wk
            tdomain1 <- tdomain
        }
        else {
            mf2 <- mf.wk
            yy2 <- yy
            mf.p2 <- mf.p.wk
            mt.p2 <- mt.p.wk
            matx.p2 <- matx.p.wk
            random2 <- random.wk
            type2 <- type.wk
            tdomain2 <- tdomain
        }
    }
    ## Generate sub-basis
    nobs <- dim(mf1)[1]
    if (is.null(id.basis)) {
        if (is.null(nbasis)) nbasis <- max(30,ceiling(10*nobs^(2/9)))
        if (nbasis>=nobs) nbasis <- nobs
        if (!is.null(seed)) set.seed(seed)
        id.basis <- sample(nobs,nbasis,prob=cnt)
    }
    else {
        if (max(id.basis)>nobs|min(id.basis)<1)
            stop("gss error in sshzd2d1: id.basis out of range")
        nbasis <- length(id.basis)
    }
    ## Generate terms
    if (symmetry) {
        nhzd <- 1
        if (dim(mf2)[2]!=dim(mf1)[2])
            stop("gss error in sshzd2d1: variables in parallel formulas must match")
        mf1.wk <- mf1
        mf2.wk <- mf2
        names(mf1.wk) <- names(mf2)
        names(mf2.wk) <- names(mf1)
        tdomain1 <- c(min(tdomain1,tdomain2),max(tdomain1,tdomain2))
        type1[[yy1$tname]][[2]] <- tdomain1
        type2[[yy2$tname]][[2]] <- tdomain1
        term1 <- mkterm(rbind(mf1,mf2.wk),type1)
        term2 <- mkterm(rbind(mf2,mf1.wk),type2)
        mf1 <- rbind(mf1,mf2.wk)
        yy1.sv <- yy1
        yy1$start <- c(yy1$start,yy2$start)
        yy1$end <- c(yy1$end,yy2$end)
        yy1$status <- c(yy1$status,yy2$status)
        id.basis.wk <- c(id.basis,id.basis+nobs)
        if (!is.null(mf.p1)) {
            if (is.null(mf.p2)||(dim(mf.p2)[2]!=dim(mf.p1)[2]))
                stop("gss error in sshzd2d1: variables in parallel formulas must match")
            matx.p1 <- rbind(matx.p1,matx.p2)
        }
        if (!is.null(random1)) {
            if (is.null(random2)||(dim(random2$z)[2]!=dim(random1$z)[2]))
                stop("gss error in sshzd2d1: variables in parallel formulas must match")
            random1$z <- rbind(random1$z,random2$z)
        }
        if (!is.null(cnt)) cnt.wk <- c(cnt,cnt)
        else cnt.wk <- NULL
    }
    else {
        nhzd <- 2
        term1 <- mkterm(mf1,type1)
        term2 <- mkterm(mf2,type2)
        id.basis.wk <- id.basis
        cnt.wk <- cnt
    }
    ## Fit marginal hazard models
    for (ii in 1:nhzd) {
        ## Extract model components
        mf <- switch(ii,mf1,mf2)
        yy <- switch(ii,yy1,yy2)
        term <- switch(ii,term1,term2)
        mf.p <- switch(ii,mf.p1,mf.p2)
        mt.p <- switch(ii,mt.p1,mt.p2)
        matx.p <- switch(ii,matx.p1,matx.p2)
        random <- switch(ii,random1,random2)
        tdomain <- switch(ii,tdomain1,tdomain2)
        if (rho=="weibull") partt <- switch(ii,part1,part2)
        ## Finalize id.basis
        nobs <- length(yy$status)
        id.basis.wk <- id.basis.wk[yy$status[id.basis.wk]]
        nbasis <- length(id.basis.wk)
        id.wk <- NULL
        nT <- sum(yy$status)
        for (i in 1:nbasis) {
            id.wk <- c(id.wk,(1:nT)[(1:nobs)[yy$status]%in%id.basis.wk[i]])
        }
        ## Generate Gauss-Legendre quadrature
        nmesh <- 200
        quad <- gauss.quad(nmesh,tdomain)
        ## set up partial terms
        if (!is.null(mf.p)) {
            for (lab in colnames(mf.p)) mf[,lab] <- mf.p[,lab]
            lab.p <- labels(mt.p)
            matx.p <- scale(matx.p)
            center.p <- attr(matx.p,"scaled:center")
            scale.p <- attr(matx.p,"scaled:scale")
            part <- list(mt=mt.p,center=center.p,scale=scale.p)
        }
        else part <- lab.p <- NULL
        ## Obtain unique covariate observations
        tname <- yy$tname
        xnames <- names(mf)
        xnames <- xnames[!xnames%in%tname]
        if (length(xnames)) {
            xx <- mf[,xnames,drop=FALSE]
            if (!is.null(part)) xx <- cbind(xx,matx.p)
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
        else stop("gss error in sshzd2d1: missing covariate")
        ## calculate rho
        if (is.null(cnt.wk)) yy$cnt <- rep(1,nobs)
        else yy$cnt <- cnt.wk
        if (rho=="marginal") {
            rho.wk <- sshzd(Surv(end,status,start)~end,data=yy,
                            id.basis=id.basis.wk,weights=cnt,alpha=2)
            rho.qd <- hzdcurve.sshzd(rho.wk,quad$pt)
            rhowk <- hzdcurve.sshzd(rho.wk,yy$end[yy$status])
        }
        if (rho=="weibull") {
            y.wk <- cbind(yy$end,yy$status,yy$start)
            form <- as.formula(paste("y.wk~",paste(xnames,collapse="+")))
            rho.wk <- gssanova(form,family="weibull",partial=partt,data=mf,
                               id.basis=id.basis.wk,weights=cnt,alpha=2)
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
            if (is.null(cnt.wk)) qd.wt[,x.ind[i]] <- qd.wt[,x.ind[i]]+wk
            else qd.wt[,x.ind[i]] <- qd.wt[,x.ind[i]]+cnt.wk[i]*wk
        }
        if (is.null(cnt.wk)) qd.wt <- quad$wt*qd.wt/nobs
        else qd.wt <- quad$wt*qd.wt/sum(cnt.wk)
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
            xy.basis <- mf[id.basis.wk,vlist]
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
        if (!is.null(part)) {
            s <- cbind(s,matx.p[yy$status,])
            int.s <- c(int.s,t(matx.p[!x.dup.ind,])%*%apply(qd.wt,2,sum))
        }
        ## generate int.z
        if (!is.null(random)) random$int.z <- t(qd.z)%*%apply(qd.wt,2,sum)
        ## Check s rank
        if (!is.null(s)) {
            nnull <- dim(s)[2]
            if (qr(s)$rank<nnull)
                stop("gss error in sshzd2d1: unpenalized terms are linearly dependent")
        }
        ## Fit the model
        Nobs <- ifelse(is.null(cnt.wk),nobs,sum(cnt.wk))
        if (!is.null(cnt.wk)) cntt <- cnt.wk[yy$status]
        else cntt <- NULL
        z <- msphzd1(s,r,id.wk,Nobs,cntt,int.s,int.r,rhowk,random,prec,maxiter,alpha,skip.iter)
        ## cfit
        if (!is.null(random)) rhowk <- rhowk*exp(-random$z%*%z$b)
        if (!is.null(cnt.wk)) cfit <- sum(cntt*rhowk)/Nobs/sum(qd.wt)
        else cfit <- sum(rhowk)/Nobs/sum(qd.wt)
        ## Brief description of model terms
        desc <- NULL
        for (label in term$labels)
            desc <- rbind(desc,as.numeric(c(term[[label]][c("nphi","nrk")])))
        if (!is.null(part)) {
            desc <- rbind(desc,matrix(c(1,0),length(lab.p),2,byrow=TRUE))
        }
        desc <- rbind(desc,apply(desc,2,sum))
        if (is.null(part)) rownames(desc) <- c(term$labels,"total")
        else rownames(desc) <- c(term$labels,lab.p,"total")
        colnames(desc) <- c("Unpenalized","Penalized")
        ## Return the results
        hzd <- c(list(call=match.call(),mf=mf,cnt=cnt.wk,terms=term,desc=desc,
                      alpha=alpha,tname=tname,xnames=xnames,tdomain=tdomain,cfit=cfit,
                      quad=quad,x.pt=x.pt,qd.wt=qd.wt,id.basis=id.basis.wk,partial=part,
                      lab.p=lab.p,random=random,skip.iter=skip.iter),z)
        hzd$se.aux$v <- sqrt(Nobs)*hzd$se.aux$v
        class(hzd) <- c("sshzd1","sshzd")
        if (ii==1) hzd1 <- hzd
        else hzd2 <- hzd
    }
    ## Finalize hzd2
    if (symmetry) {
        hzd2 <- hzd1
        hzd2$terms <- term2
        names(hzd2$mf) <- hzd2$xnames <- c(names(mf2),names(mf.p2))
        hzd2$tname <- yy2$tname
        hzd2$xnames <- names(hzd2$x.pt) <- hzd2$xnames[!hzd2$xnames%in%hzd2$tname]
        hzd2$lab.p <- labels(mt.p2)
        if (!is.null(hzd2$partial)) hzd2$partial$mt <- mt.p2
        desc <- NULL
        for (label in term2$labels)
            desc <- rbind(desc,as.numeric(c(term2[[label]][c("nphi","nrk")])))
        if (!is.null(part)) {
            desc <- rbind(desc,matrix(c(1,0),length(hzd2$lab.p),2,byrow=TRUE))
        }
        desc <- rbind(desc,apply(desc,2,sum))
        if (is.null(part)) rownames(desc) <- c(term2$labels,"total")
        else rownames(desc) <- c(term2$labels,hzd2$lab.p,"total")
        colnames(desc) <- c("Unpenalized","Penalized")
        hzd2$desc <- desc
        yy1 <- yy1.sv
    }
    ## Estimate copula
    nobs <- length(yy1$status)
    s1 <- min(hzd1$tdomain)
    s2 <- min(hzd2$tdomain)
    u1 <- u2 <- NULL
    for (i in 1:nobs) {
        u1 <- c(u1,survexp.sshzd(hzd1,yy1$end[i],data[i,,drop=FALSE],s1))
        u2 <- c(u2,survexp.sshzd(hzd2,yy2$end[i],data[i,,drop=FALSE],s2))
    }
    v1 <-v2 <- NULL
    if (any(yy1$start&yy2$start)) {
        for (i in 1:nobs) {
            v1 <- c(v1,survexp.sshzd(hzd1,yy1$start[i],data[i,,drop=FALSE],s1))
            v2 <- c(v2,survexp.sshzd(hzd2,yy2$start[i],data[i,,drop=FALSE],s2))
        }
    }
    cens <- as.numeric(!yy1$status)+as.numeric(!yy2$status)*2
    if (!is.null(v1)) trun <- cbind(v1,v2)
    else trun <- NULL
    copu <- sscopu2(cbind(u1,u2),cens,trun,symmetry,id.basis=id.basis)
    ## Return the fits
    obj <- list(call=match.call(),symmetry=symmetry,alpha=alpha,
                id.basis=id.basis,hzd1=hzd1,hzd2=hzd2,copu=copu)
    class(obj) <- "sshzd2d"
    obj
}
