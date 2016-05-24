#'multivariate Skew-Normal probability density function
#'
#'
#'@param x p x n data matrix with n the number of observations and
#'p the number of dimensions
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
#'@param Log logical flag for returning the log of the probability density
#'function. Defaults is \code{TRUE}.
#'
#'@seealso mvnpdf, mmvsnpdfC
#'
#'@importFrom stats pnorm
#'
#'@export
#'
#'@examples
#'
#'mvnpdf(x=matrix(1.96), mean=0, varcovM=diag(1), Log=FALSE)
#'dnorm(1.96)
#'mvsnpdf(x=matrix(rep(1.96,1), nrow=1, ncol=1),
#'       xi=c(0), psi=c(0), sigma=diag(1),
#'       Log=FALSE
#')
#'
#'mvsnpdf(x=matrix(rep(1.96,2), nrow=2, ncol=1),
#'       xi=c(0, 0), psi=c(1, 1), sigma=diag(2)
#')
#'
#'N=50000#00
#'Yn <- rnorm(n=N, mean=0, sd=1)
#'
#'Z <- rtruncnorm(n=N, a=0, b=Inf, mean=0, sd=1)
#'eps <- rnorm(n=N, mean=0, sd=1)
#'psi <- 10
#'Ysn <- psi*Z + eps
#'
#'nu <- 1.5
#'W <- rgamma(n=N, shape=nu/2, rate=nu/2)
#'Yst=Ysn/sqrt(W)
#'
#'library(reshape2)
#'library(ggplot2)
#'data2plot <- melt(cbind.data.frame(Ysn, Yst))
#'#pdf(file="ExSNST.pdf", height=5, width=4)
#'p <- (ggplot(data=data2plot)
#'      + geom_density(aes(x=value, fill=variable, alpha=variable), col="black")#, lwd=1.1)
#'      + theme_bw()
#'      + xlim(-15,100)
#'      + theme(legend.position="bottom")
#'      + scale_fill_manual(values=alpha(c("#F8766D", "#00B0F6"),c(0.2,0.45)),
#'                          name =" ",
#'                          labels=c("Y~SN(0,1,10)      ", "Y~ST(0,1,10,1.5)")
#'      )
#'      + scale_alpha_manual(guide=FALSE, values=c(0.25, 0.45))
#'      + xlab("Y")
#'      + ylim(0,0.08)
#'      + ylab("Density")
#'      + guides(fill = guide_legend(override.aes = list(colour = NULL)))
#'      + theme(legend.key = element_rect(colour = "black"))
#')
#'p
#'#dev.off()
#'
#'
mvsnpdf <- function(x, xi, sigma, psi, Log=TRUE){


    if(is.null(x) | is.null(xi) | is.null(sigma) | is.null(psi)){
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
        #recovering original paremters
        omega <- sigma + tcrossprod(psi)
        omegaInv <- solve(omega)
        smallomega <- diag(sqrt(diag(omega)))
        alph <- (smallomega%*%omegaInv%*%psi
                  /as.vector(sqrt(1-crossprod(psi,omegaInv)%*%psi)))

        if(dim(omega)[1]!=dim(omega)[2]){
            stop("omega is not a square matrix")
        }
        if(dim(omega)[1]!=p){
            stop("omega is of the wrong size")
        }
        part1 <- log(2) + mvnpdf(x, mean=xi, varcovM=omega, Log=TRUE)
        part2 <- stats::pnorm(t(alph)%*%diag(1/sqrt(diag(omega)))%*%(x0))
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

        part1 <- log(2) + mvnpdf(x, mean=xi, varcovM=omega, Log=TRUE)
        part2 <- mapply(FUN=function(a, o, x){
            stats::pnorm(crossprod(a,diag(1/sqrt(diag(o))))%*%(x))},
            x=x0, o=omega, a=alph)

    }

    res <- part1 + log(part2)
    if (!Log){
        res <- exp(part1)*part2
    }
    return(res)

}

