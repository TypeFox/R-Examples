lle_swissrole <-
function(N=1500,k=10,ss=FALSE,p=0.5,reg=2,iLLE=FALSE,v=0.8){

	#set dimensions
	n <- 3
	m <- 2
	
	#generate swiss roll data
	X <- 0*matrix(0,N,n)
	tt <- (3*pi/2)*(1+2*runif(N))
	height <- 21*runif(N)
	X[,1] <- tt*cos(tt) 
	X[,2] <- height 
	X[,3] <- tt*sin(tt)
	
	#perform lle
	res <- lle(X,m,k,reg=reg,ss=ss,p=p,iLLE=iLLE,v=v)
	
	Y <- res$Y
	X <- res$X
	choise <- res$choise
	
	#plot
	if( ss==0 ) col <- (tt-min(tt)) else col <- (tt[choise]-min(tt[choise])) 
	col <- col/max(col)*200
	plot_lle(Y,X,0,col)

}

