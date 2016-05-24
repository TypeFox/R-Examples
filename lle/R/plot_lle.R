plot_lle <-
function(Y,X,print=FALSE,col=3,name=as.numeric(Sys.time()),angle=60,inter=FALSE){
	
	#interactive plot
	if( inter ){
		require(rgl)
		plot3d( Y, col=col)
		return( NULL )
	}
	
	#get dimensions
	require( scatterplot3d )
	N <- dim(X)[1]
	n <- dim(X)[2]
	m <- dim(Y)[2]
	if( is.null(m) ) m <- 1
	
	#set plot parameters
	pch <- 19
	cex <- 1
	
	#set colour palette of 600 rainbow colours
	palette( rainbow(600) )
	
	#set colour vector
	if( length(col) == 1 & is.numeric( col ) ) col <- seq(0,1,length=N)*100*col
	
	
	#print direction of output file (if used)
	if( print ){
		cat( "graphics saved to",getwd(),"\n" )
		pdf( file=paste(name, ".pdf", sep="" ), width=800, height=480 )
	} else dev.new()

	#generate splitscreen
	#when outputing to file, splitscreen generates warning, which is harmless
	options("warn"=-1)
	if( n<4 ) split.screen(c(1,2))
	options("warn"=0)	
	
	#plot original data (depending on its dimension)
	if( n<4 ) screen(1)
	if( n==2 ){ 
		plot( X[,1], X[,2], col=col, main=paste("raw data n=",n,", N=",N,sep=""), pch=pch, cex=cex )
	} else if( n==3 ){
		scatterplot3d( X[,1], X[,2], X[,3], color=col, main=paste("raw data n=",n,", N=",N,sep=""), pch=pch, cex.symbols=cex, angle=angle )
	} else {
		plot( c(1,3),c(1,3), type="n" )
		text(2,2,"data with dimension n>3\n cannot be displayed", cex=1)
	}
	
	#plot embedded data (depending on its dimension)
	if( n<4 ) screen(2)
	if( m==1 ){
		plot( Y, col=col, main=paste("embedded data m=",m,", N=",N,sep=""), pch=pch, cex=cex )
	} else if( m==2 ){
		plot( Y[,2], Y[,1], col=col, main=paste("embedded data m=",m,", N=",N,sep=""), pch=pch, cex=cex )
	} else if( m==3 ){
		scatterplot3d( Y[,3], Y[,2], Y[,1], color=col, main=paste("embedded data m=",m,", N=",N,sep=""), pch=pch, cex.symbols=cex, angle=angle )
	} else {
		plot( c(1,3),c(1,3), type="n" )
		text(2,2,"data with dimension n>3\n cannot be displayed", cex=1)
	}
	
	if( print==1 ) dev.off()	
}

