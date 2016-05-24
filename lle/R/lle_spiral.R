lle_spiral <- 
function(){
	
	#set dimensions/parameters
	n <- 3
	m <- 1
	N <- 600
	k <- 5
	reg <- 2
	ss <- FALSE
	p <- 0.5
	iLLE <- FALSE
	v <- 0.8
	
	#generate elementary vectors
	v1 <- c(1,rep(0, each=n-1)) 
	v2 <- c(0,1,rep(0, each=n-2))
	v3 <- c(0,0,1,rep(0,each=n-3))
	dt <- 0.005	
	
	#prepare rotation matrices
	phix <- 0/180*pi #rotation x
	phiy <- 0/180*pi #rotation y
	phiz <- 0/180*pi #rotation z
	Rx <- matrix( c( 1, 0, 0, 0, cos(phix), -sin(phix), 0, sin(phix), cos(phix) ), 3 )
	Ry <- matrix( c( cos(phiy), 0, sin(phiy), 0, 1, 0, -sin(phiy), 0, cos(phiy) ), 3 )
	Rz <- matrix( c( cos(phiz), -sin(phiz), 0, sin(phiz), cos(phiz), 0, 0, 0, 1 ), 3 )
	
	#prepare splitscreen-plot
	plotdim <- c(3,3)
	split.screen(plotdim)
	scr <- 1
		
	for( T in seq(16,1.5,length=prod(plotdim)) ){
		
		#set parameters
		om <- 16*pi/T
		X <- 0*matrix(rnorm(n*N),nrow=N)
			
		#generate spiral data
		for( i in (1:N)){
			t=dt*i
			r=t
			X[i,]= r*v1*cos(om*t)+ r*v2*sin(om*t) + v3*t^3/30
		}
		
		#rotation
		X <- t(Rx%*%t(X))
		X <- t(Ry%*%t(X))
		X <- t(Rz%*%t(X))
	
		#perform lle
		cat( "------------\nT=",T,"\n")
		res <- lle(X,m,k,reg,ss=ss,p=p,id=TRUE,iLLE=iLLE)
		
		#plot
		if( (scr+2)%%3==0 ) zt=NULL else zt="" 
		if( scr>6 ) xt=NULL else xt=""
		if( scr%%3==0 ) yt=NULL else yt=""
		screen(scr)
		scatterplot3d(X,cex.symbols=0.5, mar=rep(.8,4), angle=10, xlab="", ylab="", zlab="", 
			type="p", lwd=2, highlight.3d=TRUE, cex.axis=0.8,
			main="", xlim=c(-4,4), ylim=c(-4,4),
			x.ticklabs=xt, y.ticklabs=yt, z.ticklabs=zt )
		if( scr < 8) text(-1.2,2,"m=1") else text(-1.2,2,"m=2") 
	
		scr <- scr + 1
	}
}

