## Calculate prediction and Bayesian SE from ssllrm objects
predict.ssllrm <- function (object,x,y=object$qd.pt,odds=NULL,se.odds=FALSE,...)
{
    if (class(object)!="ssllrm")
        stop("gss error in predict.ssllrm: not a ssllrm object")
    if ("random"%in%colnames(x)) {
        zz <- x$random
        x$random <- NULL
    }
    else zz <- NULL
    if (!all(sort(object$xnames)==sort(colnames(x))))
        stop("gss error in predict.ssllrm: mismatched x variable names")
    if (!all(sort(object$ynames)==sort(colnames(y))))
        stop("gss error in predict.ssllrm: mismatched y variable names")
    mf <- object$mf
    term <- object$term
    qd.pt <- object$qd.pt
    nmesh <- dim(qd.pt)[1]
    y.id <- NULL
    for (i in 1:dim(y)[1]) {
        if (!sum(duplicated(rbind(qd.pt,y[i,object$ynames,drop=FALSE]))))
            stop("gss error in predict.ssllrm: y value is out of range")
        wk <- FALSE
        for (j in 1:nmesh) {
            if (sum(duplicated(rbind(qd.pt[j,],y[i,object$ynames])))) y.id <- c(y.id,j)
        }
    }
    if (!is.null(odds)) {
        if (length(y.id)-length(odds))
            stop("gss error in predict.ssllrm: odds is of wrong length")
        if (!max(odds)|sum(odds))
            stop("gss error in predict.ssllrm: odds is not a contrast")
        if (sum(duplicated(y.id)))
            stop("gss error in predict.ssllrm: duplicated y in contrast")
        qd.pt <- qd.pt[y.id,,drop=FALSE]
    }
    ## Generate s, and r
    nobs <- dim(x)[1]
    nmesh <- dim(qd.pt)[1]
    nbasis <- length(object$id.basis)
    nnull <- length(object$d)
    nZ <- length(object$b)
    s <- NULL
    r <- array(0,c(nmesh,nbasis,nobs))
    nu <- nq <- 0
    for (label in term$labels) {
        vlist <- term[[label]]$vlist
        x.list <- object$xnames[object$xnames%in%vlist]
        y.list <- object$ynames[object$ynames%in%vlist]
        xy.basis <- mf[object$id.basis,vlist]
        qd.xy <- data.frame(matrix(0,nmesh,length(vlist)))
        names(qd.xy) <- vlist
        qd.xy[,y.list] <- qd.pt[,y.list]
        if (length(x.list)) xx <- x[,x.list,drop=FALSE]
        else xx <- NULL
        nphi <- term[[label]]$nphi
        nrk <- term[[label]]$nrk
        if (nphi) {
            phi <- term[[label]]$phi
            for (i in 1:nphi) {
                nu <- nu+1
                if (is.null(xx)) {
                    s.wk <- phi$fun(qd.xy[,,drop=TRUE],nu=i,env=phi$env)
                    wk <- matrix(s.wk,nmesh,nobs)
                }
                else {
                    wk <- NULL
                    for (j in 1:nobs) {
                        qd.xy[,x.list] <- xx[rep(j,nmesh),]
                        wk <- cbind(wk,phi$fun(qd.xy,i,phi$env))
                    }
                }
                s <- array(c(s,wk),c(nmesh,nobs,nu))
            }
        }
        if (nrk) {
            rk <- term[[label]]$rk
            for (i in 1:nrk) {
                nq <- nq+1
                if (is.null(xx)) {
                    r.wk <- rk$fun(qd.xy[,,drop=TRUE],xy.basis,nu=i,env=rk$env,out=TRUE)
                    r <- r + as.vector(10^object$theta[nq]*r.wk)
                }
                else {
                    wk <- NULL
                    for (j in 1:nobs) {
                        qd.xy[,x.list] <- xx[rep(j,nmesh),]
                        wk <- array(c(wk,rk$fun(qd.xy,xy.basis,i,rk$env,TRUE)),
                                    c(nmesh,nbasis,j))
                    }
                    r <- r + 10^object$theta[nq]*wk
                }
            }
        }
    }
    ## random effects
    if (nZ) {
        nz <- object$Random$sigma$env$nz
        if (is.null(zz)) z.wk <- matrix(0,nobs,nz)
        else z.wk <- as.matrix(zz)
        if (dim(z.wk)[2]!=nz)
            stop("gss error in predict.ssllrm: x$random is of wrong dimension")
        z <- nlvl <- NULL
        for (ylab in object$ynames) {
            y.wk <- mf[,ylab]
            lvl.wk <- levels(y.wk)
            nlvl.wk <- length(lvl.wk)
            nlvl <- c(nlvl,nlvl.wk)
            z.aux <- diag(1,nlvl.wk-1)
            z.aux <- rbind(z.aux,rep(-1,nlvl.wk-1))
            rownames(z.aux) <- lvl.wk
            pt.wk <- qd.pt[,ylab]
            for (i in 1:(nlvl.wk-1)) {
                for (j in 1:nmesh) {
                    z <- cbind(z,z.aux[pt.wk[j],i]*z.wk)
                }
            }
        }
        z <- aperm(array(z,c(nobs,nz,nmesh,nZ/nz)),c(3,2,4,1))
        z <- array(z,c(nmesh,nZ,nobs))
    }
    ## return
    if (is.null(odds)) {
        pdf <- NULL
        for (j in 1:nobs) {
            wk <- matrix(r[,,j],nmesh,nbasis)%*%object$c
            if (nnull) wk <- wk + matrix(s[,j,],nmesh,nnull)%*%object$d
            if (nZ) wk <- wk + matrix(z[,,j],nmesh,nZ)%*%object$b
            wk <- exp(wk)
            pdf <- cbind(pdf,wk/sum(wk))
        }
        return(t(pdf[y.id,]))
    }
    else {
        s.wk <- r.wk <- z.wk <- 0
        for (i in 1:length(odds)) {
            r.wk <- r.wk + odds[i]*r[i,,]
            if (nnull) s.wk <- s.wk + odds[i]*s[i,,]
            if (nZ) z.wk <- z.wk + odds[i]*z[i,,]
        }
        s.wk <- matrix(s.wk,nobs,nnull)
        r.wk <- t(matrix(r.wk,nbasis,nobs))
        z.wk <- t(matrix(z.wk,nZ,nobs))
        rs <- cbind(r.wk,z.wk,s.wk)
        if (!se.odds) as.vector(rs%*%c(object$c,object$b,object$d))
        else {
            fit <- as.vector(rs%*%c(object$c,object$b,object$d))
            se.fit <- .Fortran("hzdaux2",
                               as.double(object$se.aux$v), as.integer(dim(rs)[2]),
                               as.integer(object$se.aux$jpvt),
                               as.double(t(rs)), as.integer(dim(rs)[1]),
                               se=double(dim(rs)[1]), PACKAGE="gss")[["se"]]
            return(list(fit=fit,se.fit=se.fit))
        }
    }
}
