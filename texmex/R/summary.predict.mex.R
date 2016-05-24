summary.predict.mex <- function( object, mth, probs=c( .05, .5, .95 ), ... ){
    if ( is.R() ) stdev <- function( x ) sqrt( var( x ) )
    if ( missing( mth ) ) mth <- object$mth

    if (!is.null(object$replicates)){
        res <- t( sapply( object$replicates , function ( x ) apply( x, 2, mean ) ) )
    }
    else {
        res <- object$data$simulated
    }

    # Create summary function
    sumfun <- if (!is.null(object$replicates)){
                function( x , probs){
                    c( mean=mean( x ), se=stdev( x ) , quantile( x, probs=probs ) )
                }
              }
              else {
                function(x, probs){
                    c(mean=mean(x), quantile(x, probs=probs))
                }
              }

    ans <- apply( res, 2, sumfun, probs ) # Summary of expected values
    dn <- paste( dimnames( object$data$simulated)[[ 2 ]] ,"|", names( object$data$simulated )[[ 1 ]] , ">Q",100*object$pqu, sep="" )
    dimnames( ans )[[ 2 ]] <- dn

    # Get the threshold exceedence probabilities
    if (!is.null(object$replicates)){
        thres <- t(sapply(1:(dim(object$replicates[[1]])[[2]]) ,
                             function(i, x, mth){
                                 x <- sapply(x , function(x, i) x[, i], i=i)
                                 mth <- mth[i]
                                 apply(x , 2, function(x, mth)
                                                   mean(x > mth),
                                       mth = mth)
                             }, x=object$replicates, mth = mth))
        thres <- apply( thres, 1, mean)
        thres <- matrix( thres, nrow=1 )
    }
    else {
        thres <- sapply(1:dim(object$data$simulated)[[2]], function(i, x, mth){ mean(x[,i] > mth[i]) },
                            x = object$data$simulated, mth=mth)
        thres <- matrix(thres, nrow=1)
    }

    wn <- dimnames( object$data$simulated )[[ 2 ]][ 1 ]
    wth <- paste( "Q", 100*object$pqu, sep = "" )

    dn <- paste( "P(", dimnames( object$data$simulated )[[ 2 ]] , ">", signif(mth, ...),"|", wn, ">", wth, ")", sep = "" )

    dimnames( thres ) <- list( "", dn )
	
    ans <- list(ans=ans, thres=thres, call=object$call, pqu=object$pqu ,
                B = length( object$replicates ),
                which = names( object$data$simulated )[[ 1 ]],
                statistic=deparse( substitute( statistic ) )
               )

    oldClass(ans) <- "summary.predict.mex"
    ans
}

print.summary.predict.mex <- function( x, ... ){
    print( x$call, ... )

    if (x$B > 0) cat( "\nResults from", x$B, "bootstrap runs.\n" )

    cat( paste( "\nConditioned on ", x$which, " being above its ", 100*x$pqu, "th percentile.\n\n", sep = "" ) )
    cat( "\nConditional Mean and Quantiles:\n\n" )
    print( signif(x$ans,3), ... )

    cat( "\nConditional probability of threshold exceedance:\n\n" )

    print( signif(x$thres,3), ... )

    invisible()
}