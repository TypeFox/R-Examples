Fstat.fd <- function(y,yhat,argvals=NULL) {    

# observed, predicted and where to evaluate

    if( is.numeric(yhat) ){ yhat = as.vector(yhat) }

    if( (is.vector(y) & !is.vector(yhat)) | (is.fd(y) &!is.fd(yhat)) ){
        stop("y and yhat must both be either scalars or functional data objects.")
    }


    if( is.fd(y) ){
        rangeobs = y$basis$range
        rangehat = yhat$basis$range

        if( !prod(rangeobs == rangehat) ){
            stop("y and yhat do not have the same range")
        }


        if(is.null(argvals)){
            argvals = seq(rangeobs[1],rangeobs[2],length.out=101)
        }

        yvec = eval.fd(argvals,y)
        yhatvec = eval.fd(argvals,yhat)

        F = apply(yhatvec,1,var)/apply( (yvec-yhatvec)^2,1,mean)

    }
    else{
        yvec = y
        yhatvec = yhat
        F = var(yhatvec)/mean( (yvec-yhatvec)^2 )
    }

    return( list(F=F,argvals=argvals) )
}
