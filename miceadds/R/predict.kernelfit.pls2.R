
###########################################################
# prediction for kernelpls.fit2 function
predict.kernelpls.fit2 <- function( object , X , ... ){
    xmm <- matrix( object$Xmeans , nrow=nrow(X) ,ncol=ncol(X) , byrow=TRUE) 
    xc <- X - xmm
    # calculate projected scores
    scores <- xc %*% object$projection
    return(scores)
        }
###########################################################