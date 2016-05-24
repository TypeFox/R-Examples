## S3 method
project <- function (object,...) UseMethod("project")
## Calculate Kullback-Leibler projection from ssanova objects
project.ssanova <- function(object,include,...)
{
    if (class(object)[1]=="ssanova0")
        stop("gss error: square error projection is not implemented for ssanova0")
    nobs <- nrow(object$mf)
    nxi <- length(object$id.basis)
    labels.p <- object$lab.p
    ## evaluate full model
    mf <- object$mf
    yy <- predict(object,mf)
    wt <- model.weights(object$mf)
    if (!is.null(wt)) wt.wk <- sqrt(wt)
    offset <- model.offset(object$mf)
    if (!is.null(object$random)) {
        if (is.null(offset)) offset <- 0
        offset <- offset + object$random$z%*%object$b
    }
    if (!is.null(offset)) yy <- yy - offset
    ## extract terms in subspace
    s <- matrix(1,nobs,1)
    r <- NULL
    theta <- NULL
    nq.wk <- nq <- 0
    for (label in object$terms$labels) {
        if (label=="1") next
        if (label%in%labels.p) next
        x <- object$mf[,object$term[[label]]$vlist]
        x.basis <- object$mf[object$id.basis,object$term[[label]]$vlist]
        nphi <- object$term[[label]]$nphi
        nrk <- object$term[[label]]$nrk
        if (nphi) {
            phi <- object$term[[label]]$phi
            for (i in 1:nphi) {
                if (!any(label==include)) next
                s <- cbind(s,phi$fun(x,nu=i,env=phi$env))
            }
        }
        if (nrk) {
            rk <- object$term[[label]]$rk
            for (i in 1:nrk) {
                nq.wk <- nq.wk + 1
                if (!any(label==include)) next
                nq <- nq + 1
                theta <- c(theta,object$theta[nq.wk])
                r <- array(c(r,rk$fun(x,x.basis,nu=i,env=rk$env,out=TRUE)),
                           c(nobs,nxi,nq))
            }
        }
    }
    if (!is.null(object$partial)) {
        matx.p <- model.matrix(object$partial$mt,mf)[,-1,drop=FALSE]
        matx.p <- scale(matx.p)
        for (label in labels.p) {
            if (label%in%include) s <- cbind(s,matx.p[,label])
        }
    }
    ## calculate projection
    my.ls <- function(theta1=NULL) {
        if (!nq) {
            q <- matrix(0)
            sr <- cbind(s,0)
        }
        else {
            theta.wk <- 1:nq
            theta.wk[fix] <- theta[fix]
            if (nq-1) theta.wk[-fix] <- theta1
            sr <- 0
            for (i in 1:nq) sr <- sr + 10^theta.wk[i]*r[,,i]
            q <- 10^(-5)*sr[object$id.basis,]
            sr <- cbind(s,sr)
        }
        nn <- ncol(as.matrix(sr))
        nnull <- nn-nxi
        if (!is.null(wt)) {
            sr <- wt.wk*sr
            yy.wk <- wt.wk*yy
        }
        else yy.wk <- yy
        z <- .Fortran("reg",
                      as.double(sr), as.integer(nobs), as.integer(nnull),
                      as.double(q), as.integer(nxi), as.double(yy.wk),
                      as.integer(4),
                      double(1), double(1), double(1), dc=double(nn),
                      as.double(.Machine$double.eps),
                      double(nn*nn), double(nn), as.integer(rep(0,nn)),
                      double(max(nobs,nn)), integer(1), integer(1),
                      PACKAGE="gss")["dc"]
        assign("yhat",sr%*%z$dc,inherits=TRUE)
        if (!is.null(wt)) sum(wt*(yy-yhat/wt.wk)^2)/sum(wt)
        else mean((yy-yhat)^2)
    }
    cv.wk <- function(theta) cv.scale*my.ls(theta)+cv.shift
    ## initialization
    if (nq) {
        r.wk <- 0
        for (i in 1:nq) r.wk <- r.wk + 10^theta[i]*r[,,i]
        if (is.null(s)) theta.wk <- 0
        else theta.wk <- log10(sum(s^2)/ncol(s)/sum(r.wk^2)*nxi) / 2
        theta <- theta + theta.wk
        tmp <- NULL
        for (i in 1:nq) tmp <- c(tmp,10^theta[i]*sum(r[cbind(object$id.basis,1:nxi,i)]))
        fix <- rev(order(tmp))[1]
    }
    ## projection
    yhat <- NULL
    if (nq>1) {
        if (object$skip.iter) kl <- my.ls(theta[-fix])
        else {
            ## scale and shift cv
            tmp <- abs(my.ls(theta[-fix]))
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
            zz <- nlm(cv.wk,theta[-fix],stepmax=.5,ndigit=7)
            if (zz$code>3)
                warning("gss warning in project.ssanova: theta iteration fails to converge")
            kl <- my.ls(zz$est)
        }
    }
    else kl <- my.ls()
    if (!is.null(wt)) {
        yhat <- yhat/wt.wk
        ymean <- sum(wt*yy)/sum(wt)
        kl0 <- sum(wt*(yy-ymean)^2)/sum(wt)
        kl <- sum(wt*(yy-yhat)^2)/sum(wt)
        kl1 <- sum(wt*(ymean-yhat)^2)/sum(wt)
    }
    else {
        kl0 <- mean((yy-mean(yy))^2)
        kl <- mean((yy-yhat)^2)
        kl1 <- mean((mean(yy)-yhat)^2)
    }
    list(ratio=kl/kl0,kl=kl,check=(kl+kl1)/kl0)
}
