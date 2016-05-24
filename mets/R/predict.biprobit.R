##' @export
predict.biprobit <- function(object,newdata,X,Z,which=NULL,fun=NULL,type,...) {
    if (missing(newdata)) newdata <- data.frame(1)
    if (missing(X)) {
        ff <- object$formula; ff[2] <- NULL
        X <- model.matrix(ff,newdata)
    }
    if (missing(Z)) Z <- model.matrix(object$rho.formula,newdata)    
    p <- coef(object)

    h <- function(p) log(p/(1-p)) ## logit
    ih <- function(z) 1/(1+exp(-z)) ## expit
    if (!missing(type)) {
        h <- asin; ih <- sin      
        if (is.null(type)) {
            h <- ih <- identity
        }
        if (is.list(type)) {
            h <- type[[1]]; ih <- type[[2]]
        }
    }
        
    prob <- function(p,which=NULL,fun=NULL,...) {
        blen <- object$model$blen    
        beta1 <- p[seq(blen)]
        if (object$model$eqmarg) {
            beta2 <- beta1
        } else {
            beta2 <- beta1[-seq(blen/2)]
            beta1 <- beta1[seq(blen/2)]
        }    
        gamma <- p[-seq(blen)]
        m1 <- X%*%beta1
        m2 <- X%*%beta2
        r <- object$SigmaFun(gamma,cor=TRUE,Z=Z)
        pp <- data.frame(mu1=m1,mu2=m2,rho=r$rho)        
        p11 <- with(pp, pmvn(lower=c(0,0),upper=c(Inf,Inf),mu=cbind(mu1,mu2),sigma=cbind(rho),cor=TRUE))
        p10 <- with(pp, pmvn(lower=c(0,-Inf),upper=c(Inf,0),mu=cbind(mu1,mu2),sigma=cbind(rho),cor=TRUE))
        p01 <- with(pp, pmvn(lower=c(-Inf,0),upper=c(0,Inf),mu=cbind(mu1,mu2),sigma=cbind(rho),cor=TRUE))
        p00 <- with(pp, pmvn(lower=c(-Inf,-Inf),upper=c(0,0),mu=cbind(mu1,mu2),sigma=cbind(rho),cor=TRUE))
        p1 <- p11+p10
        p2 <- p11+p01
        res <- cbind(p11,p10,p01,p00,p1,p2,as.matrix(pp))
        if (!is.null(fun)) return(fun(res))
        if (!is.null(which)) {
            tr.idx <- base::which(which<7)
            res <- res[,which,drop=FALSE]
            if (length(tr.idx)>0) res[,tr.idx] <- h(res[,tr.idx])
            return(structure(res,tr.idx=tr.idx))
        }
        return(res)
    }
    pp <- prob(p,which=which,fun=fun)
    if (!is.null(fun)) {
        res <- estimate(object,prob,fun=fun)$coefmat
        if (!missing(newdata) && nrow(newdata)==nrow(res)) {
            suppressWarnings(res <- cbind(res,parameter=rownames(res),newdata))                   
        }
        return(res)
    }        
    if (!is.null(which)) {
        res <- estimate(object,prob,which=which)$coefmat[,c(1,3,4),drop=FALSE]
        tr.idx <- (attributes(pp)$tr.idx-1)*nrow(pp)
        if (length(tr.idx)>0) {
            for (ii in tr.idx) {
                idx <- ii+seq(nrow(pp))
                res[idx,] <- ih(res[idx,])
            }
        }
        rownames(res) <- rep(colnames(pp),each=nrow(pp))
        pp <- res
    }    
    if (!missing(newdata)) {
        suppressWarnings(pp <- cbind(pp,parameter=rownames(pp),newdata))
    }
    return(pp)
}
