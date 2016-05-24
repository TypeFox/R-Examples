#'multivariate skew-t  probability density function
#'
#'
#'@param x \code{p x n} data matrix with \code{n} the number of observations and
#'\code{p} the number of dimensions
#'
#'@param xi mean vector or list of mean vectors (either a vector,
#'a matrix or a list)
#'
#'@param sigma variance-covariance matrix or list of variance-covariance
#'matrices (either a matrix or a list)
#'
#'@param psi skew parameter vector or list of skew parameter vectors
#'(either a vector, a matrix or a list)
#'
#'@param df a numeric vector or a list of the degrees of freedom
#'(either a vector or a list)
#'
#'@param Log logical flag for returning the log of the probability density
#'function. Defaults is \code{TRUE}.
#'
#'@seealso mvtpdf, mvsnpdf, mmvstpdfC, mvstlikC
#'
#'@importFrom stats pt
#'
#'@export
#'
#'@examples
#'mvstpdf(x=matrix(rep(1.96,2), nrow=2, ncol=1),
#'       xi=c(0, 0), psi=c(1, 1), sigma=diag(2),
#'       df=100000000, Log=FALSE
#')
#'mvsnpdf(x=matrix(rep(1.96,2), nrow=2, ncol=1),
#'       xi=c(0, 0), psi=c(1, 1), sigma=diag(2),
#'       Log=FALSE
#')
#'mvstpdf(x=matrix(rep(1.96,2), nrow=2, ncol=1),
#'       xi=c(0, 0), psi=c(1, 1), sigma=diag(2),
#'       df=100000000
#')
#'mvsnpdf(x=matrix(rep(1.96,2), nrow=2, ncol=1),
#'       xi=c(0, 0), psi=c(1, 1), sigma=diag(2)
#')
#'
#'
mvstpdf <- function(x, xi, sigma, psi, df, Log=TRUE){

    if(is.null(x) | is.null(xi) | is.null(sigma) | is.null(psi)| is.null(df)){
        stop("some arguments are empty")
    }

    if(!is.matrix(x)){
        stop("x should be a matrix")
    }
    n <- dim(x)[2]
    p <- dim(x)[1]

    if(!is.list(xi)){
        if(is.null(xi)){
            stop("xi is empty")
        } else if(is.vector(xi) && length(xi)==p){
            x0 <- x-xi
        } else if(is.matrix(xi) && ncol(xi)==n){
            x0 <- x-xi
        } else{
            stop("wrong input for xi")
        }
    }else{
        x0 <- lapply(xi, function(v){x - v})
    }


    if(is.matrix(sigma)){
        #recovering original parameters
        omega <- sigma + tcrossprod(psi)
        omegaInv <- solve(omega)
        smallomega <- diag(sqrt(diag(omega)))
        alph <- (smallomega%*%omegaInv%*%psi
                 /as.vector(sqrt(1-crossprod(psi,omegaInv)%*%psi)))
        Qy <- apply(X=x0, MARGIN=2, FUN=function(v){crossprod(v,omegaInv)%*%v})

        if(dim(omega)[1]!=dim(omega)[2]){
            stop("omega is not a square matrix")
        }
        if(dim(omega)[1]!=p){
            stop("omega is of the wrong size")
        }
        part1 <- log(2) + mvtpdf(x, mean=xi, varcovM=omega, df=df, Log=TRUE)
        part2 <- stats::pt(q=(t(alph)%*%diag(1/sqrt(diag(omega)))%*%(x0)*
                           sqrt((df+p)/(df+Qy))),
                    df=df+p)
    }
    else{
        if(!is.list(sigma)){
            sigma <- apply(X=sigma, MARGIN=3, FUN=list)
            sigma <- lapply(sigma, FUN='[[', 1)
            x0 <- apply(X=x0, MARGIN=2, FUN=list)
            x0 <- lapply(x0, FUN='[[', 1)
            psi <- apply(X=psi, MARGIN=2, FUN=list)
            psi <- lapply(psi, FUN='[[', 1)
        }
        omega <- mapply(FUN=function(s,ps){s + tcrossprod(ps)},
                        s=sigma, ps=psi, SIMPLIFY=FALSE)
        omegaInv <- lapply(X=omega, FUN=solve)
        alph <- mapply(FUN=function(o, oI, ps){
            diag(sqrt(diag(o)))%*%oI%*%ps/sqrt(1-crossprod(ps,oI)%*%ps)[1,1]},
            o=omega, oI=omegaInv, ps=psi, SIMPLIFY=FALSE
        )
        if(is.matrix(x0[[1]])){
            Qy <- mapply(FUN=function(xx, oI){
                apply(X=xx, MARGIN=2,FUN=function(v){crossprod(v,oI)%*%v})},
                xx=x0, oI=omegaInv, SIMPLIFY=FALSE)
        }else{
            Qy <- mapply(FUN=function(v, oI){crossprod(v,oI)%*%v},
                v=x0, oI=omegaInv, SIMPLIFY=FALSE)
        }
        part1 <- log(2) + mvtpdf(x, mean=xi, varcovM=omega, df=df, Log=TRUE)
        part2 <- mapply(FUN=function(a, o, Q, d, x){
            stats::pt(q=(crossprod(a,diag(1/sqrt(diag(o))))%*%(x)*
                      sqrt((d+p)/(d+Q))),
               df=d+p)},
            x=x0, o=omega, Q=Qy, a=alph, d=df)

    }

    res <- part1 + log(part2)
    if (!Log){
        res <- exp(part1)*part2
    }
    return(res)

}

