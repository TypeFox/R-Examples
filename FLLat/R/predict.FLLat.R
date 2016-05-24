predict.FLLat <- function(object,newY=NULL,thresh=10^(-4),maxiter.T=100,...) {

    ## Error checking parameters.
    if (!inherits(object,"FLLat")) {
        stop("'object' must be of class 'FLLat'")
    }
    if (!all(length(thresh)==1,thresh>0)) {
        stop("'thresh' must be a real number > 0")
    }
    if (!all(length(maxiter.T)==1,maxiter.T>0)) {
        stop("'maxiter.T' must be an integer > 0")
    }

    ## The estimated features from object.
    B <- object$Beta

    ## If newY not given, return fitted Y values and Theta from object.
    if (is.null(newY)) {

        ## The estimated weights from object.
        pred.T <- object$Theta

        ## Fitted Y values.
        pred.Y <- B%*%pred.T
       
        return(list("pred.Y"=pred.Y,"Theta"=pred.T,"niter"=object$niter,
                    "rss"=object$rss))

    } else {

        if (!all(is.matrix(newY),is.double(newY),nrow(newY)==nrow(B))) {
            stop("'newY' must be a numeric matrix with the same number of rows as 'object$Beta'")
        }

        ## Setting weight constraint.
        sT <- 1

        ## Initializing Theta.
        old.T <- matrix(0,nrow=ncol(B),ncol=ncol(newY))

        ## Calculating predicted weights for newY.
        result <- .Call(TLatL2CR,newY,B,old.T,as.double(thresh),
                        as.integer(maxiter.T),as.double(sT))

        ## Predicted Y values.
        pred.Y <- B%*%result$Theta
        
        return(c("pred.Y"=list(pred.Y),result))

    }

}
