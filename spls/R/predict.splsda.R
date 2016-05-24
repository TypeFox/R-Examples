
# predict of either plsda or splsda

predict.splsda <-
function( object, newx, type = c("fit","coefficient"),
        fit.type = c("class","response"), ... )
{
    # newx: matrix of predictors
    # type = "fit" or "coefficient"
    
    type <- match.arg(type)  
    fit.type <- match.arg(fit.type)  
    class.fit <- object$class.fit
    classifier <- object$classifier
    W <- object$W
    A <- object$A
    T <- object$T
    x <- object$x
    p <- ncol( x )
    ngroups <- object$ngroups
    cutoff <- 0.5
    
    if ( type=="fit" )
    {        
        if ( missing(newx) )
        {
            if ( classifier=='lda' )
            {
                pred <- predict( class.fit )$class
            }
            if ( classifier=='logistic' )
            {
                ngroups <- object$ngroups
                if ( ngroups > 2 )
                {
                    pred <- as.numeric( as.vector( predict( class.fit ) ) )
                }
                if ( ngroups == 2 )
                {
                    pred <- predict( class.fit, type='response' )
                    if ( fit.type=='class' )
                    {
                        pred.class = pred > cutoff
                        pred <- as.numeric( pred.class )
                    }
                }
            }
        } else
        {   
            if ( ncol(newx)==p ) { newx <- newx[,A] }
            
            newx <- scale( newx, object$meanx[A], object$normx[A] )
            
            # latent components for test data
            
            T.test <- data.frame( newx %*% W )
            colnames(T.test) <- colnames(T)
            
            # prediction
            
            if ( classifier=='lda' )
            {
                pred <- predict( class.fit, newdata=T.test )$class
            }
            if ( classifier=='logistic' )
            {
                ngroups <- object$ngroups
                if ( ngroups > 2 )
                { pred <- predict( class.fit, newdata=T.test ) }
                if ( ngroups == 2 )
                {
                    pred <- predict( class.fit, type='response', newdata=T.test )
                    if ( fit.type=='class' )
                    {
                        pred.class = pred > cutoff
                        pred <- factor( as.numeric( pred.class ) )
                    }
                }
            }
        }
    }
    if ( type=="coefficient" )
    {
        coef.T1 <- as.matrix( coef( class.fit ) )
        
        if ( classifier=="logistic" )
        {    
            if ( ngroups == 2 )
            {
                beta0hat <- coef.T1[ 1, ]
                beta1hat <- W %*% coef.T1[ -1, ]
                betahat <- matrix( 0, (p+1), 1 )
                betahat[ 1, ] <- beta0hat
                betahat[ (A+1), ] <- beta1hat
            }
            if ( ngroups > 2 )
            {
                coef.T1 <- t( coef.T1 )
                beta0hat <- coef.T1[ 1, ]
                beta1hat <- W %*% coef.T1[ -1, ]
                betahat <- matrix( 0, (p+1), (ngroups-1) )
                betahat[ 1, ] <- beta0hat
                betahat[ (A+1), ] <- beta1hat
                
                colnames(betahat) <- paste('class 0 vs.',1:(ngroups-1))
            }
            if ( !is.null(colnames(x)) )
            {
                rownames(betahat) <- 1:nrow(betahat)
                rownames(betahat)[1] <- 'intercept'
                rownames(betahat)[2:nrow(betahat)] <- colnames(x)
            } else
            {
                rownames(betahat) <- c( 0, paste("x",1:p,sep="") )
                rownames(betahat)[1] <- 'intercept'
            }
        }
        
        if ( classifier=="lda" )
        {
            betahat <- matrix( 0, p, (ngroups-1) )
            betahat[ A, ] <- W %*% coef.T1
            
            colnames(betahat) <- colnames(coef.T1)
            
            if ( !is.null(colnames(x)) )
            {
                rownames(betahat) <- colnames(x)
            } else
            {
                rownames(betahat) <- paste("x",1:p,sep="")
            }
        }
        
        pred <- betahat
    }
    
    invisible(pred)
}

"coef.splsda" <-
function( object, ... )
{
    predict.splsda( object, type="coefficient", ... )
}
