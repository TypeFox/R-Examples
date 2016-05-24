## Evaluate hazard estimate
predict.sscox <- function (object,newdata,se.fit=FALSE,
                           include=c(object$terms$labels,object$lab.p),...)
{
    if (class(object)!="sscox") stop("gss error in predict.sscox: not a sscox object")
    nnew <- nrow(newdata)
    nbasis <- length(object$id.basis)
    nnull <- length(object$d)
    nz <- length(object$b)
    nn <- nbasis + nnull + nz
    labels.p <- object$lab.p
    ## Extract included terms
    if (!is.null(object$d)) s <- matrix(0,nnew,nnull)
    r <- matrix(0,nnew,nbasis)
    for (label in include) {
        if (label%in%labels.p) next
        xx <- object$mf[object$id.basis,object$terms[[label]]$vlist]
        xnew <- newdata[,object$terms[[label]]$vlist]
        nphi <- object$terms[[label]]$nphi
        nrk <-  object$terms[[label]]$nrk
        if (nphi) {
            iphi <- object$terms[[label]]$iphi
            phi <-  object$terms[[label]]$phi
            for (i in 1:nphi) {
                s[,iphi+(i-1)] <- phi$fun(xnew,nu=i,env=phi$env)
            }
        }
        if (nrk) {
            irk <- object$terms[[label]]$irk
            rk <- object$terms[[label]]$rk
            for (i in 1:nrk) {
                r <- r + 10^object$theta[irk+(i-1)]*
                  rk$fun(xnew,xx,nu=i,env=rk$env,out=TRUE)
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
            if (label%in%include) s[,nu] <- matx.p[,label]
        }
    }
    ## random effects
    if (nz) {
        if (is.null(newdata$random)) z.wk <- matrix(0,nnew,nz)
        else z.wk <- newdata$random
        rs <- cbind(r,z.wk,s)
    }
    else rs <- cbind(r,s)
    if (!se.fit) as.vector(exp(rs%*%c(object$c,object$b,object$d)))
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
