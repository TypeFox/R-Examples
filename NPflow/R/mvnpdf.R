#'multivariate-Normal probability density function
#'
#'
#'@param x p x n data matrix with n the number of observations and
#'p the number of dimensions
#'
#'@param mean mean vector or list of mean vectors (either a vector,
#'a matrix or a list)
#'
#'@param varcovM variance-covariance matrix or list of variance-covariance
#'matrices (either a matrix or a list)
#'
#'@param Log logical flag for returning the log of the probability density
#'function. Defaults is \code{TRUE}.
#'
#'@seealso mvnpdf, mmvnpdfC
#'
#'@author Boris P. Hejblum
#'
#'@export
#'
#'@examples
#'
#'mvnpdf(x=matrix(1.96), mean=0, varcovM=diag(1), Log=FALSE)
#'dnorm(1.96)
#'
#'mvnpdf(x=matrix(rep(1.96,2), nrow=2, ncol=1),
#'       mean=c(0, 0), varcovM=diag(2), Log=FALSE
#')
#'
mvnpdf <- function(x, mean, varcovM, Log=TRUE){
    if(!is.matrix(x)){
        stop("x should be a matrix")
    }
    n <- ncol(x)
    p <- nrow(x)


    if(!is.list(mean)){
        if(is.null(mean)){
            stop("mean is empty")
        } else if(is.vector(mean) && length(mean)==p){
            x0 <- x-mean
        } else if(is.matrix(mean) && ncol(mean)==n){
            x0 <- x-mean
        } else{
            stop("wrong input for mean")
        }
    }else{

        x0 <- lapply(mean, function(v){x - v})

    }

    if(is.null(varcovM)){
        stop("varcovM is empty")
    }


    if(is.matrix(varcovM)){
        if(dim(varcovM)[1]!=dim(varcovM)[2]){
            stop("varcovM is not a square matrix")
        }
        if(dim(varcovM)[1]!=p){
            stop("varcovM is of the wrong size")
        }

        if(is.vector(x0)){
            x0 <- matrix(x0, ncol=1)
        }

        Rinv = backsolve(chol(varcovM),x=diag(p))
        xRinv <- matrix(apply(X=x0, MARGIN=2, FUN=crossprod, y=Rinv))
        logSqrtDetvarcovM <- sum(log(diag(Rinv)))

        quadform <- apply(X=xRinv, MARGIN=2, FUN=crossprod)
        y <- (-0.5*quadform + logSqrtDetvarcovM -p*log(2*pi)/2)

        #         dMvn <- function(X,mu,Sigma) {
        #             k <- ncol(X)
        #             rooti <- backsolve(chol(Sigma),diag(k))
        #             quads <- colSums((crossprod(rooti,(t(X)-mu)))^2)
        #             return(exp(-(k/2)*log(2*pi) + sum(log(diag(rooti))) - .5*quads))
        #         }
    }else if(!is.list(mean)){
        if(!is.list(varcovM)){
            varcovM <- apply(X=varcovM, MARGIN=3, FUN=list)
            varcovM <- lapply(varcovM, FUN='[[', 1)
        }
        if(!is.list(x0)){
            x0 <- apply(X=x0, MARGIN=2, FUN=list)
            x0 <- lapply(x0, FUN='[[', 1)
        }

        likelihood <- function(x0, varcovM){
            p <- length(x0)
            Rinv = backsolve(chol(varcovM),x=diag(p))
            xRinv <- x0 %*% Rinv
            logSqrtDetvarcovM <- sum(log(diag(Rinv)))

            quadform <- tcrossprod(xRinv)
            y <- (-0.5*quadform + logSqrtDetvarcovM -p*log(2*pi)/2)
        }
        y <-mapply(FUN=likelihood, x0, varcovM)


    }else{
        R <- lapply(varcovM, FUN=chol)
        Rinv <- lapply(X=R, FUN=backsolve, x=diag(p))
        xRinv <- mapply(FUN=function(x,y){apply(X=x,MARGIN=2, FUN=crossprod, y=y)}, x=x0, y=Rinv, SIMPLIFY=FALSE)
        logSqrtDetvarcovM <- lapply(X=Rinv, FUN=function(X){sum(log(diag(X)))})
        quadform <- lapply(X=xRinv, FUN=function(x){apply(X=x, MARGIN=2, FUN=crossprod)})
        y <- mapply(FUN=function(x,y){(-0.5*x + y -p*log(2*pi)/2)},
                    x=quadform, y=logSqrtDetvarcovM)
    }

    if(!Log){
        y <- exp(y)
    }

    return(y)

}