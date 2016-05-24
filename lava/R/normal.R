##' @export
rmvn <- function(n,mu=rep(0,ncol(sigma)),sigma=diag(nrow=length(mu))*.5+.5,...) {
    PP <- with(svd(sigma), v%*%diag(sqrt(d),ncol=length(d))%*%t(u))
    res <- matrix(rnorm(ncol(sigma)*n),ncol=ncol(sigma))%*%PP
    if (NROW(mu)==nrow(res) && NCOL(mu)==ncol(res)) return(res+mu)
    return(res+cbind(rep(1,n))%*%mu)
}

##' @export
dmvn <- function(x,mu,sigma,log=FALSE,nan.zero=TRUE,norm=TRUE,...) {
    if (length(sigma)==1) {
        k <- 1
        isigma <- structure(cbind(1/sigma),det=as.vector(sigma))

    } else {
        k <- ncol(sigma)
        isigma <- Inverse(sigma)
    }
    if (!missing(mu)) {
        if (NROW(mu)==NROW(x) && NCOL(mu)==NCOL(x)) {
            x <- x-mu
        } else {
            x <- t(t(x)-mu)
        }
    }
    logval <- -0.5*(base::log(2*base::pi)*k+
                    base::log(attributes(isigma)$det)+
                    rowSums((x%*%isigma)*x))
    if (nan.zero) logval[is.nan(logval)] <- -Inf
    if (log) return(logval)
    return(exp(logval))
}


normal_method.lvm <- "nlminb2"

normal_objective.lvm <- function(x,p,data,weight2=NULL,indiv=FALSE,...) {
    if (!exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) runif(1)
    save.seed <- get(".Random.seed", envir = .GlobalEnv)
    on.exit(assign(".Random.seed", save.seed, envir = .GlobalEnv))
    set.seed(1)
    y.idx <- lava::index(x)$endo.idx
    y <- lava::endogenous(x)
    ord <- lava::ordinal(x)
    status <- rep(0,length(y))
    if (exists("binary.lvm")) {
        bin <- match(do.call("binary.lvm",list(x=x)),y)
        if (length(bin)>0) status[bin] <- 2
    }
    status[match(ord,y)] <- 2

    Table <- length(y)==length(ord)
    if (Table) {
        pat <- mets::fast.pattern(data,categories=max(data)+1)
        data <- pat$pattern
        colnames(data) <- y
    }

    mu <- predict(x,data=data,p=p)
    S <- attributes(mu)$cond.var
    class(mu) <- "matrix"
    thres <- matrix(0,nrow=length(y),max(1,attributes(ord)$K-1)); rownames(thres) <- y
    for (i in seq_len(length(attributes(ord)$fix))) {
        nn <- names(attributes(ord)$idx)[i]
        ii <- attributes(ord)$idx[[nn]]
        val <- (attributes(mu)$e[ii])
        thres[nn,seq_len(length(val))] <-
            cumsum(c(val[1],exp(val[-1])))
    }

    yl <- yu <- as.matrix(data[,y,drop=FALSE])
    if (!inherits(yl[1,1],c("numeric","integer","logical")) ||
        !inherits(yu[1,1],c("numeric","integer","logical")))
        stop("Unexpected data (normal_objective)")
    if (!is.null(weight2)) {
        yu[,colnames(weight2)] <- weight2
        status[match(colnames(weight2),y)] <- 1
    }
    l <- mets::loglikMVN(yl,yu,status,mu,S,thres)

    if (Table) {
        l <- l[pat$group+1]
        ##data <- data[pat$group+1,]
        ##l <- l+runif(length(l),0,0.001)
    }

    if (indiv) return(-l)
    return(-sum(l))
}

normal_logLik.lvm <- function(object,p,data,weight2=NULL,...) {
    res <- -normal_objective.lvm(x=object,p=p,data=data,weight2=weight2,...)
    return(res)
}

normal_gradient.lvm <- function(x,p,data,weight2=NULL,indiv=FALSE,...) {
    if  (is.null(ordinal(x))) {
        D <- deriv.lvm(x,p=p)
        M <- moments(x,p)
        Y <- as.matrix(data[,manifest(x)])
        mu <- t(M$xi)%x%rep(1,nrow(Y))
        ss <- -mets::scoreMVN(Y,mu,M$C,D$dxi,D$dS)
        if (!indiv) return(colSums(ss))
        return(ss)
    }
    if (indiv) {
        return(numDeriv::jacobian(function(p0) normal_objective.lvm(x,p=p0,data=data,weight2=weight2,indiv=TRUE,...),p,method=lava.options()$Dmethod))
    }
    numDeriv::grad(function(p0) normal_objective.lvm(x,p=p0,data=data,weight2=weight2,...),p,method=lava.options()$Dmethod)
}

normal_hessian.lvm <- function(x,p,n,...) {
    ##return(numDeriv::jacobian(function(p0) normal_gradient.lvm(x,p=p0,data=data,indiv=FALSE,...),p,method=lava.options()$Dmethod))
    dots <- list(...); dots$weight <- NULL
    do.call("information", c(list(x=x,p=p,n=n),dots))
    ## S <- normal_gradient.lvm(x,p=p,data=data,indiv=TRUE,...)
    ## J <- t(S)%*%S
    ## attributes(J)$grad <- colSums(S)
    ## return(J)
}
