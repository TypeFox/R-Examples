## Calculate prediction and Bayesian SE from ssanova0 objects
predict.ssanova0 <- function(object,newdata,se.fit=FALSE,
                             include=c(object$terms$labels,object$lab.p),...)
{
    nnew <- dim(newdata)[1]
    nobs <- length(object$c)
    nnull <- length(object$d)
    labels.p <- object$lab.p
    ## Extract included terms
    term <- object$terms
    philist <- rklist <- NULL
    s <- q <- NULL
    nq <- 0
    for (label in include) {
        if (label=="1") {
            philist <- c(philist,term[[label]]$iphi)
            s <- cbind(s,rep(1,len=nnew))
            next
        }
        if (label%in%labels.p) next
        if (label=="offset") next
        xnew <- newdata[,term[[label]]$vlist]
        x <- object$mf[,term[[label]]$vlist]
        nphi <- term[[label]]$nphi
        nrk <- term[[label]]$nrk
        if (nphi) {
            iphi <- term[[label]]$iphi
            phi <- term[[label]]$phi
            for (i in 1:nphi) {
                philist <- c(philist,iphi+(i-1))
                s <- cbind(s,phi$fun(xnew,nu=i,env=phi$env))
            }
        }
        if (nrk) {
            irk <- term[[label]]$irk
            rk <- term[[label]]$rk
            for (i in 1:nrk) {
                rklist <- c(rklist,irk+(i-1))
                nq <- nq+1
                q <- array(c(q,rk$fun(xnew,x,nu=i,env=rk$env,out=TRUE)),c(nnew,nobs,nq))
            }
        }
    }
    if (!is.null(object$partial)) {
        vars.p <- as.character(attr(object$partial$mt,"variables"))[-1]
        facs.p <- attr(object$partial$mt,"factors")
        vlist <- vars.p[as.logical(apply(facs.p,1,sum))]
        for (lab in labels.p) {
            if (lab%in%include) {
                vlist.wk <- vars.p[as.logical(facs.p[,lab])]
                vlist <- vlist[!(vlist%in%vlist.wk)]
            }
        }
        if (length(vlist)) {
            for (lab in vlist) newdata[[lab]] <- 0
        }
        matx.p <- model.matrix(object$partial$mt,newdata)[,-1,drop=FALSE]
        matx.p <- sweep(matx.p,2,object$partial$center)
        matx.p <- sweep(matx.p,2,object$partial$scale,"/")
        nu <- nnull-dim(matx.p)[2]
        for (label in labels.p) {
            nu <- nu+1
            if (label%in%include) {
                philist <- c(philist,nu)
                s <- cbind(s,matx.p[,label])
            }
        }
    }
    qq <- matrix(0,nnew,nobs)
    nq <- 0
    for (i in rklist) {
        nq <- nq + 1
        qq <- qq + 10^object$theta[i]*q[,,nq]
    }
    if (!is.null(object$w)) w <- object$w
    else w <- model.weights(object$mf)
    if (!is.null(w)) qq <- t(sqrt(w)*t(qq))
    ## Compute posterior mean
    nphi <- length(philist)
    pmean <- as.vector(qq%*%object$c)
    if (nphi) pmean <- pmean + as.vector(s%*%object$d[philist])
    if (any(include=="offset")) {
        if (is.null(model.offset(object$mf)))
            stop("gss error: no offset in the fit")
        offset <- newdata$offset
        if (is.null(offset)) offset <- newdata$"(offset)"
        if (is.null(offset)) stop("gss error: missing offset")
        pmean <- pmean + offset
    }
    if (se.fit) {
        b <- object$varht/10^object$nlambda
        ## Get cr, dr, and sms
        crdr <- getcrdr(object,t(qq))
        cr <- crdr$cr
        dr <- crdr$dr[philist,,drop=FALSE]
        sms <- getsms(object)[philist,philist]
        ## Compute posterior variance
        r <- 0
        for (label in include) {
            if (label=="1") next
            if (label%in%labels.p) next
            xnew <- newdata[,term[[label]]$vlist]
            nrk <- term[[label]]$nrk
            if (nrk) {
                irk <- term[[label]]$irk
                rk <- term[[label]]$rk
                for (i in 1:nrk) {
                    ind <- irk+(i-1)
                    r <- r + 10^object$theta[ind]*rk$fun(xnew,xnew,nu=i,env=rk$env)
                }
            }
        }
        fn2 <- function(x,n) x[1:n]%*%x[n+(1:n)]
        pvar <- r - apply(rbind(t(qq),cr),2,fn2,nobs)
        if (nphi) {
            fn1 <- function(x,sms) t(x)%*%sms%*%x
            pvar <- pvar + apply(s,1,fn1,sms)
            pvar <- pvar - 2*apply(rbind(t(s),dr),2,fn2,nphi)
        }
        pse <- as.numeric(sqrt(b*pvar))
        list(fit=pmean,se.fit=pse)
    }
    else pmean
}
