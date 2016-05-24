
# Heatmap-type MSPE plot for s*K

"cv.sgpls.multi" <-
function( x, y, fold=10, K, eta, scale.x=TRUE, plot.it=TRUE,
    br=TRUE, ftype='iden', n.core=8 )
{
    result.mat <- c()

    # data partition

    foldi <- cv.split( y, fold )

    # initialization

    x <- as.matrix(x)
    n <- nrow(x)
    p <- ncol(x)
    ip <- c(1:p)
    y <- as.matrix(y)
    q <- ncol(y)

    # eta & K pair

    eta.K.pair <- cbind( rep(eta,each=length(K)), rep(K,length(eta)) )
    eta.K.list <- split( eta.K.pair, c(1:nrow(eta.K.pair)) )

    # fit SPLSDA for given eta & K

    .fit.sgpls.multi <- function( eta.val, K.val ) {
	mspemati <- rep( 0, fold )
	Ai <- rep( 0, fold )

	for ( k in 1:fold )
	{
		# fold

		#print( paste('fold ',k,sep='') )

		omit <- foldi[[k]]
		train.x <- x[-omit,]
		train.y <- y[-omit,]
		test.x <- x[omit,]
		test.y <- y[omit,]

		sgpls.fit <- sgpls.multi( train.x, train.y, K=K.val, eta=eta.val,
			scale.x=scale.x, br=br, ftype=ftype )
		pred <- as.numeric( as.vector( predict(sgpls.fit, newx=test.x) ) )
		mspemati[ k ] <- mean( as.numeric( pred != test.y ) )
		Ai[ k ] <- mean( length(sgpls.fit$A.all) )
	}

	mspe.ij <- c( mean(mspemati), mean(Ai), eta.val, K.val )
	return(mspe.ij)
    }

    # CV MSPE estimation
    if (.Platform$OS.type == "unix") {
        result.list <- parallel::mclapply( eta.K.list,
            function(x) .fit.sgpls.multi( eta.val=x[1], K.val=x[2] ),
            mc.cores = n.core )
    } else {
        # otherwise, use usual "lapply"

        result.list <- lapply( eta.K.list,
            function(x) .fit.sgpls.multi( eta.val=x[1], K.val=x[2] ) )
    }

    result.mat <- c()
    for ( i in 1:length(result.list) ) {
    	result.mat <- rbind( result.mat, result.list[[i]] )
    }

    mspemat <- matrix( result.mat[,1], length(K), length(eta) )
    mspemat <- t(mspemat)
    rownames(mspemat) <- eta
    colnames(mspemat) <- K

    # find optimal eta & K

    cands <- result.mat[ result.mat[,1]==min(result.mat[,1]),,drop=FALSE ]
    cands <- cands[ cands[,2]==min(cands[,2]),,drop=FALSE ]
    cands <- cands[ cands[,4]==min(cands[,4]),,drop=FALSE ]
    cands <- cands[ cands[,3]==max(cands[,3]),,drop=FALSE ]

    K.opt <- cands[ , 4 ]
    eta.opt <- cands[ , 3 ]

    cat( paste('\nOptimal parameters: eta = ',eta.opt,', ',sep='') )
    cat( paste('K = ',K.opt,'\n',sep='') )

    # plot heatmap & return values

    if ( plot.it )
    { heatmap.spls( t(mspemat), xlab='K', ylab='eta', main='CV MSPE Plot', coln=16, as='n' ) }
    rownames(mspemat) <- paste('eta=',eta)
    colnames(mspemat) <- paste('K =',K)

    cv <- list( mspemat=mspemat, eta.opt=eta.opt, K.opt=K.opt )
    invisible(cv)
}
