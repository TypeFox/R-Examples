LTT.plot <- function (trees,width=1,precalc=0,bound=10^(-12),timemax=10,nmax=10,avg=FALSE) {
	if (precalc==1){ltt<-trees} else {ltt<-LTT.plot.gen(trees,bound)}
	lttavg <- ltt[[1]]	
	plot(c(0,-timemax),c(1,nmax),type="l",xlab = "Time before present", ylab = "Number of species",log="y", xaxt="n",yaxt="n",lwd=width,col="white")  	        
	axis(1,at=5*(0:-100),labels=5*(0:100))
	for (i in 2:length(ltt)){
		ltt2 <- ltt[[(i)]]
		lines(ltt2)
	}
	axis(2,at=c(1,2,5,10,50,100,500,1000,5000),labels=c(1,2,5,10,50,100,500,1000,5000),xpd=TRUE)
	if (avg==TRUE){ lines(lttavg[,1],lttavg[,2],col="red")}
	}