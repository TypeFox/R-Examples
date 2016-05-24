spcr <- function(x, y, k, lambda.B, lambda.gamma, w=0.1, xi=0.01, adaptive=FALSE, center=TRUE, scale=FALSE){
	if( !is.matrix(x) ) stop("x must be a matrix.")
	if( mode(x)!="numeric" ) stop("x must be numeric.")
	if ( !is.vector(y) ) stop("y must be a vector.")
	if( mode(y)!="numeric" ) stop("y must be numeric.")
	
#	ini.lambda.B <- ini.lambda.gamma <- ini.lambda( x=x, y=y, k=k, w=w, xi=xi )
#	if( ini.lambda.B < lambda.B ) stop("lambda.B is large. Set smaller lambda.B.")
#	if( ini.lambda.gamma < lambda.gamma ) stop("lambda.gamma is large. Set smaller lambda.gamma.")
	
	if( center==TRUE ) x <- sweep(x, 2, apply(x,2,mean))
	if( scale==TRUE ) x <- scale(x)
	
	A <- as.matrix(eigen(var(x))$vectors[ ,1:k])
	gamma0 <- mean(y)
	gamma <- rep(0, k)
	Beta <- matrix( 0, nrow(A), k )

	if( adaptive==FALSE ){
		spcr.object <- .Call( "spcr", x, y, A, Beta, gamma, gamma0, lambda.B, lambda.gamma, xi, w )
		ans <- list( loadings.B=spcr.object[[1]], gamma=spcr.object[[2]], gamma0=spcr.object[[3]], loadings.A=spcr.object[[4]], call=match.call() )
		class(ans) <- "spcr"
		ans
	} else {
		spcr.object <- .Call( "spcr", x, y, A, Beta, gamma, gamma0, lambda.B, lambda.gamma, xi, w )
		Beta <- spcr.object[[1]]
		gamma <- spcr.object[[2]]
		gamma0 <- spcr.object[[3]]
		A <- spcr.object[[4]]
		BetaWeight <- Beta/sum(abs(Beta))		
		adaspcr.object <- .Call( "adaspcr", x, y, A, Beta, gamma, gamma0, lambda.B, lambda.gamma, xi, w , BetaWeight)
		ans <- list( loadings.B=adaspcr.object[[1]], gamma=adaspcr.object[[2]], gamma0=adaspcr.object[[3]], loadings.A=adaspcr.object[[4]], call=match.call() )
		class(ans) <- "spcr"
		ans		
	}
}
