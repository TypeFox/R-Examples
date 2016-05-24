lle_sound <-
function(t=500,dt=20,k=25,reg=2,ss=FALSE,p=0.5,id=TRUE){

	#generate data from wavefile
	#result: dataset with 'dt' dimension and approx '(f/dt-t/dt)' samples
	i <- 1 
	invisible( data( lle_wave ) )
	f <- length(lle_wave)
	X <- t( matrix( lle_wave[1:t] ) )
	while( i+t+dt < f ){
		i <- i + dt
		X <- rbind( X, lle_wave[i:(i+t-1)] )
	}
	
	#get/set parameters
	n <- dim(X)[2]
	N <- dim(X)[1]
	m <- 8
	
	#perform lle
	res <- lle(X=X,m=m,k=k,reg=reg,ss=ss,p=p,id=id)
	
	Y <- res$Y
	X <- res$X
	choise <- res$choise
	
	#plot
	plot( res$id, type="l", main="found intrinsic dimension")
	lines( smooth.spline( res$id, df=5), col="red", lwd=2)
	
}

