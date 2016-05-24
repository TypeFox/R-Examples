hzdrate.sshzd <- ## Evaluate hazard estimate
function (object,x,se=FALSE,include=c(object$terms$labels,object$lab.p)) {
    if (!any(class(object)=="sshzd"))
        stop("gss error in hzdrate.sshzd: not a sshzd object")
    if (dim(object$mf)[2]==1&is.vector(x)) {
        x <- data.frame(x)
        colnames(x) <- colnames(object$mf)
    }
    if (!is.null(object$d)) s <- matrix(0,dim(x)[1],length(object$d))
    r <- matrix(0,dim(x)[1],length(object$id.basis))
    for (label in include) {
        if (label=="1") {
            iphi <- object$terms[[label]]$iphi
            s[,iphi] <- rep(1,dim(x)[1])
            next
        }
        if (label%in%object$lab.p) next
        xx <- object$mf[object$id.basis,object$terms[[label]]$vlist]
        x.new <- x[,object$terms[[label]]$vlist]
        nphi <- object$terms[[label]]$nphi
        nrk <-  object$terms[[label]]$nrk
        if (nphi) {
            iphi <- object$terms[[label]]$iphi
            phi <-  object$terms[[label]]$phi
            for (i in 1:nphi) {
                s[,iphi+(i-1)] <- phi$fun(x.new,nu=i,env=phi$env)
            }
        }
        if (nrk) {
            irk <- object$terms[[label]]$irk
            rk <- object$terms[[label]]$rk
            for (i in 1:nrk) {
                r <- r + 10^object$theta[irk+(i-1)]*
                  rk$fun(x.new,xx,nu=i,env=rk$env,out=TRUE)
            }
        }
    }
    if (!is.null(object$partial)) {
        vars.p <- as.character(attr(object$partial$mt,"variables"))[-1]
        facs.p <- attr(object$partial$mt,"factors")
        vlist <- vars.p[as.logical(apply(facs.p,1,sum))]
        for (lab in object$lab.p) {
            if (lab%in%include) {
                vlist.wk <- vars.p[as.logical(facs.p[,lab])]
                vlist <- vlist[!(vlist%in%vlist.wk)]
            }
        }
        if (length(vlist)) {
            for (lab in vlist) x[[lab]] <- 0
        }
        matx.p <- model.matrix(object$partial$mt,x)[,-1,drop=FALSE]
        matx.p <- sweep(matx.p,2,object$partial$center)
        matx.p <- sweep(matx.p,2,object$partial$scale,"/")
        nu <- length(object$d)-dim(matx.p)[2]
        for (label in object$lab.p) {
            nu <- nu+1
            if (label%in%include) s[,nu] <- matx.p[,label]
        }
    }
    if (is.null(object$random)) rs <- cbind(r,s)
    else {
        nz <- length(object$b)
        rs <- cbind(r,matrix(0,dim(x)[1],nz),s)
    }
    if (!se) as.vector(exp(rs%*%c(object$c,object$b,object$d)))
    else {
        fit <- as.vector(exp(rs%*%c(object$c,object$b,object$d)))
        se.fit <- .Fortran("hzdaux2",
                           as.double(object$se.aux$v), as.integer(dim(rs)[2]),
                           as.integer(object$se.aux$jpvt),
                           as.double(t(rs)), as.integer(dim(rs)[1]),
                           se=double(dim(rs)[1]), PACKAGE="gss")[["se"]]
        list(fit=fit,se.fit=se.fit)
    }
}

hzdcurve.sshzd <- ## Evaluate hazard curve for plotting
function (object,time,covariates=NULL,se=FALSE) {
    tname <- object$tname
    xnames <- object$xnames
    if (!any(class(object)=="sshzd"))
        stop("gss error in hzdcurve.sshzd: not a sshzd object")
    if (length(xnames)&&(!all(xnames%in%names(covariates))))
        stop("gss error in hzdcurve.sshzd: missing covariates")
    mn <- min(object$tdomain)
    mx <- max(object$tdomain)
    if ((min(time)<mn)|(max(time)>mx))
        stop("gss error in hzdcurve.sshzd: time range beyond the domain")
    if (length(xnames)) {
        xx <- covariates[,xnames,drop=FALSE]
        xy <- data.frame(matrix(0,length(time),length(xnames)+1))
        names(xy) <- c(tname,xnames)
        xy[,tname] <- time
    }
    else xx <- NULL
    if (!se) {
        if (is.null(xx))
            zz <- hzdrate.sshzd(object,time)
        else {
            zz <- NULL
            for (i in 1:dim(xx)[1]) {
                xy[,xnames] <- xx[rep(i,length(time)),]
                zz <- cbind(zz,hzdrate.sshzd(object,xy))
            }
            zz <- zz[,,drop=TRUE]
        }
        zz
    }
    else {
        if (is.null(xx))
            zz <- hzdrate.sshzd(object,time,TRUE)
        else {
            fit <- se.fit <- NULL
            for (i in 1:dim(xx)[1]) {
                xy[,xnames] <- xx[rep(i,length(time)),]
                wk <- hzdrate.sshzd(object,xy,TRUE)
                fit <- cbind(fit,wk$fit)
                se.fit <- cbind(se.fit,wk$se.fit)
            }
            zz <- list(fit=fit[,,drop=TRUE],se.fit=se.fit[,,drop=TRUE])
        }
        zz
    }
}

survexp.sshzd <- ## Compute expected survival
function(object,time,covariates=NULL,start=0) {
    tname <- object$tname
    xnames <- object$xnames
    ## Check inputs
    if (!any(class(object)=="sshzd"))
        stop("gss error in survexp.sshzd: not a sshzd object")
    if (length(xnames)&&(!all(xnames%in%names(covariates))))
        stop("gss error in survexp.sshzd: missing covariates")
    lmt <- cbind(start,time)
    if (any(lmt[,1]>lmt[,2]))
        stop("gss error in survexp.sshzd: start after follow-up time")
    nt <- dim(lmt)[1]
    if (is.null(covariates)) ncov <- 1
    else ncov <- dim(covariates)[1]
    mn <- min(object$tdomain)
    mx <- max(object$tdomain)
    if ((min(start)<mn)|(max(time)>mx))
        stop("gss error in survexp.sshzd: time range beyond the domain")
    ## Calculate
    qd.hize <- 200
    qd <- gauss.quad(2*qd.hize,c(mn,mx))
    gap <- diff(qd$pt)
    g.wk <- gap[qd.hize]/2
    for (i in 1:(qd.hize-2)) g.wk <- c(g.wk,gap[qd.hize+i]-g.wk[i])
    g.wk <- 2*g.wk
    g.wk <- c(g.wk,(mx-mn)/2-sum(g.wk))
    gap[qd.hize:1] <- gap[qd.hize+(1:qd.hize)] <- g.wk
    brk <- cumsum(c(mn,gap))
    if (is.null(covariates)) {
        zz <- NULL
        d.qd <- hzdcurve.sshzd(object,qd$pt)
        for (i in 1:nt) {
            ind <- (1:(2*qd.hize))[(qd$pt<lmt[i,2])&(qd$pt>lmt[i,1])]
            if (length(ind)) {
                wk <- sum(d.qd[ind]*qd$wt[ind])
                id.mx <- max(ind)
                if (lmt[i,2]<brk[id.mx+1])
                    wk <- wk-d.qd[id.mx]*qd$wt[id.mx]*(brk[id.mx+1]-lmt[i,2])/gap[id.mx]
                else wk <- wk+d.qd[id.mx+1]*qd$wt[id.mx+1]*(lmt[i,2]-brk[id.mx+1])/gap[id.mx+1]
                id.mn <- min(ind)
                if (lmt[i,1]<brk[id.mn])
                    wk <- wk+d.qd[id.mn-1]*qd$wt[id.mn-1]*(brk[id.mn]-lmt[i,1])/gap[id.mn-1]
                else wk <- wk-d.qd[id.mn]*qd$wt[id.mn]*(lmt[i,1]-brk[id.mn])/gap[id.mn]
            }
            else {
                if (lmt[i,1]<=qd$pt[1])
                    wk <- d.qd[1]*qd$wt[1]*(lmt[i,2]-lmt[i,1])/gap[1]
                if (lmt[i,1]>=qd$pt[2*qd.hize])
                    wk <- d.qd[2*qd.hize]*qd$wt[2*qd.hize]*(lmt[i,2]-lmt[i,1])/gap[2*qd.hize]
                if ((lmt[i,1]>qd$pt[1])&(lmt[i,1]<qd$pt[2*qd.hize])) {
                    i.wk <- min((1:(2*qd.hize))[qd$pt>lmt[i,1]])
                    if (brk[i.wk]<=lmt[i,1])
                        wk <- d.qd[i.wk]*qd$wt[i.wk]*(lmt[i,2]-lmt[i,1])/gap[i.wk]
                    if (brk[i.wk]>=lmt[i,2])
                        wk <- d.qd[i.wk-1]*qd$wt[i.wk-1]*(lmt[i,2]-lmt[i,1])/gap[i.wk-1]
                    if ((brk[i.wk]<lmt[i,2])&(brk[i.wk]>lmt[i,1]))
                        wk <- d.qd[i.wk]*qd$wt[i.wk]*(lmt[i,2]-brk[i.wk])/gap[i.wk]+
                          d.qd[i.wk-1]*qd$wt[i.wk-1]*(brk[i.wk]-lmt[i,1])/gap[i.wk-1]
                }
            }
            zz <- c(zz,wk)
        }
    }
    else {
        zz <- NULL
        for (j in 1:ncov) {
            zz.wk <- NULL
            d.qd <- hzdcurve.sshzd(object,qd$pt,covariates[j,,drop=FALSE])
            for (i in 1:nt) {
                ind <- (1:(2*qd.hize))[(qd$pt<lmt[i,2])&(qd$pt>lmt[i,1])]
                if (length(ind)) {
                    wk <- sum(d.qd[ind]*qd$wt[ind])
                    id.mx <- max(ind)
                    if (lmt[i,2]<=brk[id.mx+1])
                        wk <- wk-d.qd[id.mx]*qd$wt[id.mx]*(brk[id.mx+1]-lmt[i,2])/gap[id.mx]
                    else wk <- wk+d.qd[id.mx+1]*qd$wt[id.mx+1]*(lmt[i,2]-brk[id.mx+1])/gap[id.mx+1]
                    id.mn <- min(ind)
                    if (lmt[i,1]<brk[id.mn])
                        wk <- wk+d.qd[id.mn-1]*qd$wt[id.mn-1]*(brk[id.mn]-lmt[i,1])/gap[id.mn-1]
                    else wk <- wk-d.qd[id.mn]*qd$wt[id.mn]*(lmt[i,1]-brk[id.mn])/gap[id.mn]
                }
                else {
                    if (lmt[i,1]<=qd$pt[1])
                        wk <- d.qd[1]*qd$wt[1]*(lmt[i,2]-lmt[i,1])/gap[1]
                    if (lmt[i,1]>=qd$pt[2*qd.hize])
                        wk <- d.qd[2*qd.hize]*qd$wt[2*qd.hize]*(lmt[i,2]-lmt[i,1])/gap[2*qd.hize]
                    if ((lmt[i,1]>qd$pt[1])&(lmt[i,1]<qd$pt[2*qd.hize])) {
                        i.wk <- min((1:(2*qd.hize))[qd$pt>lmt[i,1]])
                        if (brk[i.wk]<=lmt[i,1])
                            wk <- d.qd[i.wk]*qd$wt[i.wk]*(lmt[i,2]-lmt[i,1])/gap[i.wk]
                        if (brk[i.wk]>=lmt[i,2])
                            wk <- d.qd[i.wk-1]*qd$wt[i.wk-1]*(lmt[i,2]-lmt[i,1])/gap[i.wk-1]
                        if ((brk[i.wk]<lmt[i,2])&(brk[i.wk]>lmt[i,1]))
                            wk <- d.qd[i.wk]*qd$wt[i.wk]*(lmt[i,2]-brk[i.wk])/gap[i.wk]+
                              d.qd[i.wk-1]*qd$wt[i.wk-1]*(brk[i.wk]-lmt[i,1])/gap[i.wk-1]
                    }
                }
                zz.wk <- c(zz.wk,wk)
            }
            zz <- cbind(zz,as.vector(zz.wk))
        }
        if (ncov==1) zz <- as.vector(zz)
    }
    exp(-zz)
}
