lle_scurve <-
function(N=800,k=12,ss=FALSE,p=0.5,reg=2,iLLE=FALSE,v=0.8){
	
	#set dimensions
	n <- 3
	m <- 2
	
	#generate S-curve with noise
	X <- 0*matrix(rnorm(n*N),nrow=N)
	angle <- pi*(1.5*runif(N/2)-1)
	height <- 5*runif(N)
	sd <- 0.03
	X[,1] <- c(cos(angle),-cos(angle)) + rnorm(1,0,sd)
	X[,2] <- height + rnorm(1,0,sd)
	X[,3] <- c( sin(angle),2-sin(angle) ) + rnorm(1,0,sd)
	
	#lle
	res <- lle(X,m,k,reg=reg,ss=ss,p=p,iLLE=iLLE)
	Y <- res$Y
	X <- res$X
	choise <- res$choise
	
	#plot
	col <- c(angle,angle)
	if( ss==1 ) col <- col[choise]
	col <- ((col-mean(col))/min(col)+1)*80
	plot_lle(Y,X,0,col,"",10)
}

