
# main SPLS fit

"spls" <-
function( x, y, K, eta, kappa=0.5, select="pls2", fit="simpls",
        scale.x=TRUE, scale.y=FALSE, eps=1e-4, maxstep=100, trace=FALSE )
{
    # always required to input: x, y, K, plsmethod
    # x: matrix
    # y: matrix or vector
    # K: number of latent components
    # eta: penalty parameter; 0<eta<1
    # kappa: if y multivariate; 0<kappa<=0.5
    # select: "simpls" (update X) or "pls2" (update Y)
    # fit: "simpls" or "kernelpls" or "widekernelpls" or "oscorespls"
    # scale.x: scale x?
    # scale.y: scale y?
    # eps: if y multivariate; criterion to determine whether a=c
    # maxstep: if y multivariate; force finish each iteration

    # check x, y matrix/vector

    # initialization

    x <- as.matrix(x)
    n <- nrow(x)
    p <- ncol(x)
    ip <- c(1:p)
    y <- as.matrix(y)
    q <- ncol(y)
    one <- matrix(1,1,n)

    # center & scale y & x

    mu <- one %*% y / n
    y <- scale( y, drop(mu), FALSE )
    meanx <- drop( one %*% x ) / n
    x <- scale( x, meanx, FALSE )

    if ( scale.x )
    {
        normx <- sqrt( drop( one %*% (x^2) ) / (n-1) )
        if ( any( normx < .Machine$double.eps ) )
        { stop("Some of the columns of the predictor matrix have zero variance.") }
        x <- scale( x, FALSE, normx )
    } else { normx <- rep( 1, p ) }

    if ( scale.y )
    {
        normy <- sqrt( drop( one %*% (y^2) ) / (n-1) )
        if ( any( normy < .Machine$double.eps ) )
        { stop("Some of the columns of the response matrix have zero variance.") }
        y <- scale( y, FALSE, normy )
    } else { normy <- rep( 1, q ) }

    # initilize objects

    betahat <- matrix( 0, p, q )
    betamat <- list()
    x1 <- x
    y1 <- y

    type <- correctp( x, y, eta, K, kappa, select, fit )
    eta <- type$eta
    K <- type$K
    kappa <- type$kappa
    select <- type$select
    fit <- type$fit

    # main iteration

    if ( is.null(colnames(x)) )
    { xnames <- c(1:p) } else { xnames <- colnames(x) }

    new2As <- list()
    if ( trace )
    { cat("The variables that join the set of selected variables at each step:\n") }

    for (k in 1:K)
    {
        # define Z

        Z <- t(x1) %*% y1

        # fit direction vector

        what <- spls.dv( Z, eta, kappa, eps, maxstep )

        # construct A

        A <- unique( ip[ what!=0 | betahat[,1]!=0 ] )
        new2A <- ip[ what!=0 & betahat[,1]==0 ]

        # fit pls with predictors in A

        xA <- x[,A,drop=FALSE]
        plsfit <- pls::plsr( y~xA, ncomp=min(k,length(A)),
                            method=fit, scale=FALSE )

        # update

        betahat <- matrix( 0, p, q )
        betahat[A,] <- matrix( coef(plsfit), length(A), q )
        betamat[[k]] <- betahat # for cv.spls
        pj <- plsfit$projection

        if ( select=="pls2" )
        {
            y1 <- y - x %*% betahat
        }
        if ( select=="simpls" )
        {
            pw <- pj %*% solve( t(pj) %*% pj ) %*% t(pj)
            x1 <- x
            x1[,A] <- x[,A,drop=FALSE] - x[,A,drop=FALSE] %*% pw
        }

        # print out variables that join the active set

        new2As[[k]] <- new2A

        if ( trace )
        {
            if ( length(new2A)<=10 )
            {
                cat( paste("- ",k,"th step (K=",k,"):\n",sep="") )
                cat( xnames[new2A] )
                cat( "\n" )
            } else
            {
                cat( paste("- ",k,"th step (K=",k,"):\n",sep="") )
                nlines <- ceiling(length(new2A)/10)
                for ( i in 0:(nlines-2) )
                {
                    cat( xnames[new2A[(10*i+1):(10*(i+1))]] )
                    cat( "\n" )
                }
                cat( xnames[new2A[(10*(nlines-1)+1):length(new2A)]] )
                cat( "\n" )
            }
        }
    }

    # return objects

    if ( !is.null(colnames(x)) ) { rownames(betahat) <- colnames(x) }
    if ( q>1 & !is.null(colnames(y)) ) { colnames(betahat) <- colnames(y) }

    object <- list( x=x, y=y, betahat=betahat, A=A, betamat=betamat, new2As=new2As,
                    mu=mu, meanx=meanx, normx=normx, normy=normy,
                    eta=eta, K=K, kappa=kappa, select=select, fit=fit,
                    projection=pj )
    class(object) <- "spls"
    object
}
