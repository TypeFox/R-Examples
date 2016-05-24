calc_k <- 
function(X,m,kmin=1,kmax=20,plotres=TRUE,parallel=FALSE,cpus=2,iLLE=FALSE){
	
	#set parameters
	N <- dim(X)[1]
	if( kmax>=N ) kmax <- N - 1 #more neighbourse than points doesnt make sense
	if( .Platform$OS.type=="windows" ) dev <- "nul" else dev <- "/dev/null"	
	
	#set up parallel computation
	if( parallel==TRUE) sfInit( parallel=TRUE, cpus=cpus ) else sfInit( parallel=FALSE ) 
	#require lle for every node
	options("warn"=-1)
	sfLibrary( lle )
	options("warn"=0)
	
	
	perform_calc <- function( k, X, m, iLLE=FALSE ){

		N <- dim(X)[1]

		#perform LLE
		sink( dev ) #surpress output
		Y <- lle(X,m,k,2,0,iLLE=iLLE)$Y
		sink()

		#distance matrix of original data
		Dx <- as.matrix(dist(X))

		#distance matrix of embedded data
		Dy <- as.matrix(dist(Y))

		#calculate correlation between original and embedded data for every data point
		rho <- c()
		for( i in 1:N ) rho <- c( rho, cor(Dx[i,],Dy[i,]) ) 

		#rho cant be cumulated, rho^2 can
		return( mean(1-rho^2) )
	}

	rho <- invisible( sfLapply( kmin:kmax, perform_calc, X, m, iLLE ) )
	rho <- unclass( unlist( rho ) )
	sfStop()	

	res <- data.frame( k=c(kmin:kmax), rho=rho )

	if( plotres ){
		par( mar=c(5,5,4,2)+0.1 )
		plot( res$k, res$rho, type="b", xlab="k", ylab=expression(1-rho^2), main="" )
		abline(h=min(res$rho,na.rm=TRUE),col="red")
		grid()
	} else cat( "best k:",head(res$k[order(res$rho)],3), "\n\n" )	
	return ( res )
}