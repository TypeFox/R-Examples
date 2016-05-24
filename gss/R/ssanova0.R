## Fit ssanova model
ssanova0 <- function(formula,type=NULL,data=list(),weights,subset,
                     offset,na.action=na.omit,partial=NULL,
                     method="v",varht=1,prec=1e-7,maxiter=30)
{
    ## Obtain model frame and model terms
    mf <- match.call()
    mf$type <- mf$method <- mf$varht <- mf$partial <- NULL
    mf$prec <- mf$maxiter <- NULL
    mf[[1]] <- as.name("model.frame")
    mf <- eval(mf,parent.frame())
    ## Generate terms
    term <- mkterm(mf,type)
    ## Generate s, q, and y
    nobs <- dim(mf)[1]
    s <- q <- NULL
    nq <- 0
    for (label in term$labels) {
        if (label=="1") {
            s <- cbind(s,rep(1,len=nobs))
            next
        }
        x <- mf[,term[[label]]$vlist]
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
                q <- array(c(q,rk$fun(x,x,nu=i,env=rk$env,out=TRUE)),c(nobs,nobs,nq))
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
            stop("gss error in ssanova: partial data are of wrong size")
        matx.p <- scale(matx.p)
        center.p <- attr(matx.p,"scaled:center")
        scale.p <- attr(matx.p,"scaled:scale")
        s <- cbind(s,matx.p)
        part <- list(mt=mt.p,center=center.p,scale=scale.p)
    }
    else part <- lab.p <- NULL
    ## Prepare the data
    y <- model.response(mf,"numeric")
    w <- model.weights(mf)
    offset <- model.offset(mf)
    if (!is.null(offset)) {
        term$labels <- c(term$labels,"offset")
        term$offset <- list(nphi=0,nrk=0)
        y <- y - offset
    }
    if (!is.null(w)) {
        w <- sqrt(w)
        y <- w*y
        s <- w*s
        for (i in 1:nq) q[,,i] <- w*t(w*q[,,i])
    }
    if (qr(s)$rank<dim(s)[2])
        stop("gss error in ssanova0: unpenalized terms are linearly dependent")
    if (!nq) stop("gss error in ssanova0: use lm for models with only unpenalized terms")
    ## Fit the model
    if (nq==1) {
        q <- q[,,1]
        z <- sspreg0(s,q,y,method,varht)
    }
    else z <- mspreg0(s,q,y,method,varht,prec,maxiter)
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
    obj <- c(list(call=match.call(),mf=mf,terms=term,partial=part,lab.p=lab.p,
                  desc=desc),z)
    class(obj) <- c("ssanova0","ssanova")
    obj
}

## Fit Single Smoothing Parameter REGression
sspreg0 <- function(s,q,y,method="v",varht=1)
{
    ## Check inputs
    if (is.vector(s)) s <- as.matrix(s)
    if (!(is.matrix(s)&is.matrix(q)&is.vector(y)&is.character(method))) {
        stop("gss error in sspreg: inputs are of wrong types")
    }
    nobs <- length(y)
    nnull <- dim(s)[2]
    if (!((dim(s)[1]==nobs)&(dim(q)[1]==nobs)&(dim(q)[2]==nobs)
          &(nobs>=nnull)&(nnull>0))) {
        stop("gss error in sspreg: inputs have wrong dimensions")
    }
    ## Set method for smoothing parameter selection
    code <- (1:3)[c("v","m","u")==method]
    if (!length(code)) {
        stop("gss error: unsupported method for smoothing parameter selection")
    }
    ## Call RKPACK driver DSIDR
    z <- .Fortran("dsidr0",
                  as.integer(code),
                  swk=as.double(s), as.integer(nobs),
                  as.integer(nobs), as.integer(nnull),
                  as.double(y),
                  qwk=as.double(q), as.integer(nobs),
                  as.double(0), as.integer(0), double(2),
                  nlambda=double(1), score=double(1), varht=as.double(varht),
                  c=double(nobs), d=double(nnull),
                  qraux=double(nnull), jpvt=integer(nnull),
                  double(3*nobs),
                  info=integer(1),PACKAGE="gss")
    ## Check info for error
    if (info<-z$info) {               
        if (info>0)
            stop("gss error in sspreg: matrix s is rank deficient")
        if (info==-2)
            stop("gss error in sspreg: matrix q is indefinite")
        if (info==-1)
            stop("gss error in sspreg: input data have wrong dimensions")
        if (info==-3)
            stop("gss error in sspreg: unknown method for smoothing parameter selection.")
    }
    ## Return the fit
    c(list(method=method,theta=0),
      z[c("c","d","nlambda","score","varht","swk","qraux","jpvt","qwk")])
}

## Fit Multiple Smoothing Parameter REGression
mspreg0 <- function(s,q,y,method="v",varht=1,prec=1e-7,maxiter=30)
{
    ## Check inputs
    if (is.vector(s)) s <- as.matrix(s)
    if (!(is.matrix(s)&is.array(q)&(length(dim(q))==3)
          &is.vector(y)&is.character(method))) {
        stop("gss error in mspreg: inputs are of wrong types")
    }
    nobs <- length(y)
    nnull <- dim(s)[2]
    nq <- dim(q)[3]
    if (!((dim(s)[1]==nobs)&(dim(q)[1]==nobs)&(dim(q)[2]==nobs)
          &(nobs>=nnull)&(nnull>0)&(nq>1))) {
        stop("gss error in mspreg: inputs have wrong dimensions")
    }
    ## Set method for smoothing parameter selection
    code <- (1:3)[c("v","m","u")==method]
    if (!length(code)) {
        stop("gss error: unsupported method for smoothing parameter selection")
    }
    ## Call RKPACK driver DMUDR
    z <- .Fortran("dmudr0",
                  as.integer(code),
                  as.double(s),         # s
                  as.integer(nobs), as.integer(nobs), as.integer(nnull),
                  as.double(q),         # q
                  as.integer(nobs), as.integer(nobs), as.integer(nq),
                  as.double(y),         # y
                  as.double(0), as.integer(0),
                  as.double(prec), as.integer(maxiter),
                  theta=double(nq), nlambda=double(1),
                  score=double(1), varht=as.double(varht),
                  c=double(nobs), d=double(nnull),
                  double(nobs*nobs*(nq+2)),
                  info=integer(1),PACKAGE="gss")[c("theta","info")]
    ## Check info for error
    if (info<-z$info) {               
        if (info>0)
            stop("gss error in mspreg: matrix s is rank deficient")
        if (info==-2)
            stop("gss error in mspreg: matrix q is indefinite")
        if (info==-1)
            stop("gss error in mspreg: input data have wrong dimensions")
        if (info==-3)
            stop("gss error in mspreg: unknown method for smoothing parameter selection.")
        if (info==-4)
            stop("gss error in mspreg: iteration fails to converge, try to increase maxiter")
        if (info==-5)
            stop("gss error in mspreg: iteration fails to find a reasonable descent direction")
    }
    qwk <- 10^z$theta[1]*q[,,1]
    for (i in 2:nq) qwk <- qwk + 10^z$theta[i]*q[,,i]
    ## Call RKPACK driver DSIDR
    zz <- .Fortran("dsidr0",
                   as.integer(code),
                   swk=as.double(s), as.integer(nobs),
                   as.integer(nobs), as.integer(nnull),
                   as.double(y),
                   qwk=as.double(qwk), as.integer(nobs),
                   as.double(0), as.integer(0), double(2),
                   nlambda=double(1), score=double(1), varht=as.double(varht),
                   c=double(nobs), d=double(nnull),
                   qraux=double(nnull), jpvt=integer(nnull),
                   double(3*nobs),
                   info=integer(1),PACKAGE="gss")
    ## Return the fit
    c(list(method=method,theta=z$theta),
      zz[c("c","d","nlambda","score","varht","swk","qraux","jpvt","qwk")])
}

## Obtain c & d for new y's
getcrdr <- function(obj,r)
{
    ## Check inputs
    if (is.vector(r)) r <- as.matrix(r)
    if (!(any(class(obj)=="ssanova0")&is.matrix(r))) {
        stop("gss error in getcrdr: inputs are of wrong types")
    }
    nobs <- length(obj$c)
    nnull <- length(obj$d)
    nr <- dim(r)[2]
    if (!((dim(r)[1]==nobs)&(nr>0))) {
        stop("gss error in getcrdr: inputs have wrong dimensions")
    }
    ## Call RKPACK ulitity DCRDR
    z <- .Fortran("dcrdr",
                  as.double(obj$swk), as.integer(nobs),
                  as.integer(nobs), as.integer(nnull),
                  as.double(obj$qraux), as.integer(obj$jpvt),
                  as.double(obj$qwk), as.integer(nobs),
                  as.double(obj$nlambda),
                  as.double(r), as.integer(nobs), as.integer(nr),
                  cr=double(nobs*nr), as.integer(nobs),
                  dr=double(nnull*nr), as.integer(nnull),
                  double(2*nobs), integer(1),PACKAGE="gss")[c("cr","dr")]
    ## Return cr and dr
    z$cr <- matrix(z$cr,nobs,nr)
    z$dr <- matrix(z$dr,nnull,nr)
    z
}

## Obtain var-cov matrix for unpenalized terms
getsms <- function(obj)
{
    ## Check input
    if (!any(class(obj)=="ssanova0")) {
        stop("gss error in getsms: inputs are of wrong types")
    }
    nobs <- length(obj$c)
    nnull <- length(obj$d)
    ## Call RKPACK ulitity DSMS
    z <- .Fortran("dsms",
                  as.double(obj$swk), as.integer(nobs),
                  as.integer(nobs), as.integer(nnull),
                  as.integer(obj$jpvt),
                  as.double(obj$qwk), as.integer(nobs),
                  as.double(obj$nlambda),
                  sms=double(nnull*nnull), as.integer(nnull),
                  double(2*nobs), integer(1),PACKAGE="gss")["sms"]
    ## Return the nnull-by-nnull matrix
    matrix(z$sms,nnull,nnull)
}
