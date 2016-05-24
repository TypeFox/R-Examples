
test.gMCP.correlation <- function() {
	Gm <- matrix(0,nr=4,nc=4) 
	Gm[1,3] <- 1 
	Gm[2,4] <- 1 
	Gm[3,2] <- 1 
	Gm[4,1] <- 1 
	w <- c(1/2,1/2,0,0) 
	G <- matrix2graph(Gm,w) 
	Cm <- matrix(NA,nr=4,nc=4) 
	diag(Cm) <- 1 
	Cm[1,2] <- 1/2 
	Cm[2,1] <- 1/2 
	Cm[3,4] <- 1/2 
	Cm[4,3] <- 1/2 
	p <- c(0.0131,0.1,0.012,0.01) 
	x <- unname(gMCP(G,p,corr=Cm,test="Bretz2011",alpha=0.025)@rejected)
	checkEquals(c(TRUE, FALSE, TRUE, FALSE), x)
}

test.gMCP.correlation.Alpha.Simulation <- function() {
	if (!gMCP:::tests("extended")) {
		cat("Skipping alpha level simulation.\n")
		return()
	}
	N <- 100000
	sigma <- matrix(0.9, nrow=3, ncol=3)
	diag(sigma) <- 1
	w <- c(0.4,0.4,0.2)
	G <- matrix(1, nrow=3, ncol=3)
	diag(G) <- 0
	graph <- matrix2graph(G,w)
	i <- 0
	for(k in 1:N){
		# simulate data with n obs and m-endpts			
		y=rmvnorm(1,mean=rep(0,3),sigma=sigma)
		p <- as.numeric(pnorm(y))
		res=gMCP(graph,p,corr=sigma,alpha=0.05,test="Bretz2011")
		if(min(res@adjPValues)<=0.05){i=i+1}
	}
	eps <- 0.1
	cat("Family-wise error is ",i/N, ".\n", sep="")
	checkTrue(0.05-eps<i/N && i/N<0.05+eps)
}