
# classification

"splsda" <-
function( x, y, K, eta, kappa=0.5,
        classifier=c('lda','logistic'), scale.x=TRUE, ... )
{
    y <- factor(y)
    groups <- levels(y)
    ngroups <- length(groups)
    n <- length(y)
    p <- ncol(x)
    classifier <- match.arg(classifier)

    # coding of response

    if ( ngroups > 2 )
    {
        y.code <- +1 * model.matrix(~factor(y)-1)
        colnames(y.code) <- groups
    }
    if ( ngroups == 2 )
    {
        y.code <- matrix(0,n,1)
        y.code[ y==levels(y)[2] ] <- +1
    }

    # center & scale y & x

    one <- matrix(1,1,n)

    meanx <- drop( one %*% x ) / n
    x <- scale( x, meanx, FALSE )

    if ( scale.x )
    {
        normx <- sqrt( drop( one %*% (x^2) ) / (n-1) )
        if ( any( normx < .Machine$double.eps ) )
        { stop("Some of the columns of the predictor matrix have zero variance.") }
        x <- scale( x, FALSE, normx )
    } else { normx <- rep( 1, p ) }

    # fit

    spls.fit <- spls( x, y.code, eta=eta, K=K, kappa=kappa,
                        select="simpls", scale.x=FALSE, scale.y=FALSE,
                        trace=FALSE, ... )
    A <- spls.fit$A

    # latent components

    W <- spls.fit$projection
    T <- x[,A] %*% W
    colnames(T) <- paste( 'T', c(1:K), sep='' )
    T.train <- data.frame( y, T )

    # classification

    if ( classifier=='lda' )
    { class.fit <- MASS::lda( y ~ ., data = T.train, ... ) }
    if ( classifier=='logistic' )
    {
        if ( ngroups > 2 )
        { class.fit <- nnet::multinom( y ~ ., data = T.train, trace=FALSE, ... ) }
        if ( ngroups == 2 )
        { class.fit <- glm( y ~ ., data = T.train, family=binomial, ... ) }
    }

    # return

    object <- list( spls.fit=spls.fit, class.fit=class.fit,
                    eta=eta, K=K, kappa=kappa,
                    T=T, W=W, A=A, x=x, y=y,
                    meanx=meanx, normx=normx,
                    groups=groups, ngroups=ngroups, classifier=classifier )
    class(object) <- "splsda"
    object
}
