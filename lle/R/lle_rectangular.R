lle_rectangular <-
function(N=40,k=5,v=0.9){

	#dimension
	n <- 500 
	t <- seq( 0, 1, length=n ) 
	
	#function to generate rectangular signal
	rects <- function(t, width, distance ){
		data <- t*0
		t0 <- t[1]
		t1 <- head(t,1)
		data[t>t0+distance & t<t0+distance+width] <- rnorm(1,1,0.01)
		return( data )
	}
	
	#set parameters
	width_range <- seq( 0.02, 0.35, length=N) #max(dist)+max(width) < 1 (see documentation)
	dist_range <- seq( 0.05, 0.6, length=N )
	
	#generate data
	X <- matrix( 0, ncol=n, nrow=N*N )
	for( i in 1:N ){
		for( j in 1:N ){
			w <- width_range[i]
			d <- dist_range[j]
			X[(i-1)*N+j,] <- rects( t, w, d )
		}
	}
	
	
	#perform lle
	m <- 2
	ss <- FALSE
	p <- 0.8
	reg <- 2
	id <- TRUE
	iLLE <- FALSE
	
	res <- lle(X=X,m=m,k=k,reg=reg,ss=ss,p=p,id=id,iLLE=iLLE,v=v)
	
	Y <- res$Y
	X <- res$X
	choise <- res$choise
	
	#plot
	plot_lle(Y,X,0,col=3,inter=TRUE)
}

