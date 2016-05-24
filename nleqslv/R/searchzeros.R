searchZeros <- function(x, fn, digits=4L, ... )
{
    if( !is.numeric(x) ) stop('argument x should be numeric')
    if( !is.matrix(x)  ) stop('argument x must be a matrix')
    if( any(is.na(x))  ) stop("argument x may not contain NA")
    if( !is.numeric(digits) ) stop('argument digits should be numeric')
    if( is.na(digits)  ) digits <- 4L

    N <- nrow(x)
    if( N < 1 ) stop("Matrix 'x' must have at least 1 row")

    tcode <- numeric(N)
    fnorm <- numeric(N)
    xmat  <- matrix(NA, nrow=N, ncol=ncol(x))

    kerr <- numeric(N)
    kptr <- 0

    # for k==1 check that all arguments are correct --> no try

    for ( k in seq_len(N) ){
        if( k == 1 ) {
            z <- nleqslv(x[k, ], fn=fn, ...)
        } else {
            z <- try(nleqslv(x[k, ], fn=fn, ...), silent=TRUE)
            if( inherits(z, "try-error") ) {
                kptr <- kptr+1
                kerr[kptr] <- k
                next
            }
        }
        tcode[k]  <- z$termcd
        fnorm[k]  <- norm(z$fvec,"2")^2/2 # criterion for global methods
        xmat[k, ] <- z$x
    }

    # locate index of converged trials and store corresponding starting values
    # return NULL if no full convergence obtained
    if(!any(tcode==1)) return(NULL)
    idxcvg <- which(tcode==1)
    xstartcvg <- x[idxcvg,,drop=FALSE]
    # rounded solutions for locating duplicates and remove duplicates
    xsol <- round(xmat[idxcvg,,drop=FALSE], digits)
    notdups <- !duplicated(xsol)
    xsol <- xsol[notdups,,drop=FALSE]
    solstart <- xstartcvg[notdups,,drop=FALSE]
    if( !is.null(colnames(x)) ) {
        colnames(xmat) <- colnames(x)
        colnames(xsol) <- colnames(x)
        colnames(solstart) <- colnames(x)
    }

    # order the rounded solution
    if( nrow(xsol) > 1 ) {
        zidxo <- do.call(order,split(xsol,col(xsol)))
    } else {
        zidxo <- 1
    }

    idxfatal <- if(kptr) kerr[1:kptr] else integer(0)
    idxxtol   <- which(tcode==2)
    idxnocvg <- which(tcode>2)
    # original unrounded solutions with duplicates (above) removed
    xsol <- xmat[idxcvg,,drop=FALSE][notdups,,drop=FALSE]

    # return full precision solutions ordered with rounded ordering
    res <- list(x=xsol[zidxo,,drop=FALSE], xfnorm=fnorm[idxcvg][notdups][zidxo],
                fnorm=fnorm[idxcvg], idxcvg=idxcvg, idxxtol=idxxtol,
                idxnocvg=idxnocvg, idxfatal=idxfatal,
                xstart=solstart[zidxo,,drop=FALSE],cvgstart=xstartcvg)
    res
}
