
# predict of sgpls

predict.sgpls <-
function( object, newx, type = c("fit","coefficient"),
        fit.type = c("class","response"), ... )
{
    # newx: matrix of predictors
    # type = "fit" or "coefficient"
    
    type <- match.arg(type)    
    fit.type <- match.arg(fit.type)  
    betahat <- object$betahat
    cutoff <- 0.5
    
    if ( type=="fit" )
    {        
        if ( missing(newx) )
        {
            x0 <- object$x0
        } else
        {            
            x0 <- scale( newx, object$mu, object$sigma )
        }
        
        y <- object$y
        nclass <- length(unique(y))
        
        if ( nclass == 2 )
        {
            # prediction
            
            h <- function(x) { ifelse(is.infinite(exp(x)),1,exp(x)/(1+exp(x))) }
                    
            eta.hat <- cbind( rep( 1, nrow(x0) ), x0 ) %*% betahat
            pred.prob = h( eta.hat )
            if ( fit.type=='class' )
            {
                pred.class = pred.prob > cutoff
                pred <- as.numeric( as.vector( pred.class ) )
            }
            if ( fit.type=='response' )
            { pred <- pred.prob }
        }
        if ( nclass > 2 )
        {
            # prediction
                    
            eta.hat <- cbind( rep( 1, nrow(x0) ), x0 ) %*% betahat
            mu.hat <- exp( eta.hat )
            mu.hat.sum <- 1 + apply( mu.hat, 1, sum )
            mu.hat <- mu.hat / mu.hat.sum
            mu.hat <- cbind( (1-apply(mu.hat,1,sum)), mu.hat )
            if ( fit.type=='class' )
            {
                pred.class <- apply( mu.hat, 1, which.max ) - 1
                pred <- as.numeric( as.vector( pred.class ) )
            }
            if ( fit.type=='response' )
            { pred <- mu.hat }
        }
    }
    if ( type=="coefficient" ) { pred <- betahat }
    
    invisible(pred)
}

"coef.sgpls" <-
function( object, ... )
{
    predict.sgpls( object, type="coefficient", ... )
}
