`example_model2` <- function(p, fig=FALSE){
	if(length(p)!=4){
		print("Expecting 4 Parameters")
		return(NaN)
	}
	t <- array(dim=c(4,200))
	t[1,] <- seq(0,1,length.out=200)*p[1]
	t[2,] <- seq(1,0,length.out=200)*p[2]
	t[3,] <- sin(2*pi/10*(1:200))*p[3]
	t[4,] <- rep(c(0,1,0.5,1,0),each=40)*p[4]
	if(fig){
		plot(c(0,200),c(0,max(t)),type="n")
		for(i in 1:4)
			lines(t[i,],col=i)
		legend(c(paste("P",1:4,sep=""),"sum"),lty=1,col=1:5, x=150, y=0.9*max(t))
	}
	sums<-colSums(t)
	if(fig)
		lines(sums,col=5)
	return(sums)
}

