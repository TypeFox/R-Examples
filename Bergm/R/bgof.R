bgof <- function(x,
                 directed=FALSE,
                 sample.size=100,
                 aux.iters=10000,
                 n.deg=NULL,
                 n.dist=NULL,
                 n.esp=NULL,
                 n.ideg=NULL,
                 n.odeg=NULL,
                 ...){
					
if(x$nchains > 1) x$Theta <- apply(x$Theta,2,cbind)

F <- as.matrix(x$Theta[sample(dim(x$Theta)[1],sample.size),])

if(directed==FALSE){ # undirected
	for(i in 1:sample.size){
		a <- gof.formula(x$formula,
		         coef=F[i,],
		         verbose=FALSE,
		         control=control.gof.formula(nsim=1, MCMC.burnin=aux.iters))
		if(i==1) A<-as.vector(a$pobs.deg)
		A<-cbind(A,as.vector(a$psim.deg))
		if(i==1) B<-as.vector(a$pobs.dist) 
		B<-cbind(B,as.vector(a$psim.dist))
		if(i==1) C<-as.vector(a$pobs.espart)
		C<-cbind(C,as.vector(a$psim.espart))
	}
	if(is.null(n.deg)) n.deg <- dim(A)[1]
	if(is.null(n.dist)) n.dist <- dim(B)[1]-1
	if(is.null(n.esp)) n.esp <- dim(C)[1]
	
	a5 <- apply(A[1:n.deg,-1],
	            1,quantile,probs=0.05)
	b5 <- apply(B[-(n.dist:(dim(B)[1]-1)),-1],
	            1,quantile,probs=0.05)
	c5 <- apply(C[1:n.esp,-1],
	            1,quantile,probs=0.05)
	a95 <- apply(A[1:n.deg,-1],
	            1,quantile,probs=0.95)
	b95 <- apply(B[-(n.dist:(dim(B)[1]-1)),-1],
	            1,quantile,probs=0.95)
	c95 <- apply(C[1:n.esp,-1],
	            1,quantile,probs=0.95)	
	par(mfrow=c(1,3),oma=c(0,0,3,0),mar=c(4,3,1.5,1))
	
	boxplot(as.data.frame(t(A[1:n.deg,-1])),
	        xaxt="n",
	        xlab="degree",
	        ylab="proportion of nodes")
	axis(1,seq(1,n.deg),seq(0,n.deg-1))
	lines(A[1:n.deg,1],lwd=2,col=2)
	lines(a5,col="darkgray")
	lines(a95,col="darkgray")
	
	title("Bayesian goodness-of-fit diagnostics",outer=TRUE)

	boxplot(as.data.frame(t(B[-(n.dist:(dim(B)[1]-1)),-1])),
	        xaxt="n",
	        xlab="minimum geodesic distance",
	        ylab="proportion of dyads")
	axis(1,seq(1,n.dist),labels=c(seq(1,(n.dist-1)),"NR"))
	lines(B[-(n.dist:(dim(B)[1]-1)),1],lwd=2,col=2)
	lines(b5,col="darkgray")
	lines(b95,col="darkgray")
	
	boxplot(as.data.frame(t(C[1:n.esp,-1])),
	        xaxt="n",
	        xlab="edge-wise shared partners",
	        ylab="proportion of edges")
	axis(1,seq(1,n.esp),seq(0,n.esp-1))
	lines(C[1:n.esp,1],lwd=2,col=2)
	lines(c5,col="darkgray")
	lines(c95,col="darkgray")

	out = list(sim.degree=A[,-1],
	           sim.dist=B[,-1],
	           sim.esp=C[,-1],
	           obs.degree=A[,1],
	           obs.dist=B[,1],
	           obs.esp=C[,1],
	           fun=F)
	
}else{ # directed
	
	for(i in 1:sample.size){
		a <- gof.formula(x$formula,
		         coef=F[i,],
		         verbose=FALSE,
		         GOF=~idegree+odegree+espartners+distance,
		         control=control.gof.formula(nsim=1, MCMC.burnin=aux.iters))
		if(i==1) A<-as.vector(a$pobs.ideg)
		A<-cbind(A,as.vector(a$psim.ideg))
		if(i==1) AA<-as.vector(a$pobs.odeg)
		AA<-cbind(AA,as.vector(a$psim.odeg))
		if(i==1) B<-as.vector(a$pobs.dist) 
		B<-cbind(B,as.vector(a$psim.dist))
		if(i==1) C<-as.vector(a$pobs.espart)
		C<-cbind(C,as.vector(a$psim.espart))
	}
	if(is.null(n.ideg)) n.ideg <- dim(A)[1]
	if(is.null(n.odeg)) n.odeg <- dim(AA)[1]
	if(is.null(n.dist)) n.dist <- dim(B)[1]-1
	if(is.null(n.esp)) n.esp <- dim(C)[1]
	
	a5 <- apply(A[1:n.ideg,-1],1,quantile,probs=0.05)
	aa5 <- apply(AA[1:n.odeg,-1],1,quantile,probs=0.05)
	b5 <- apply(B[-(n.dist:(dim(B)[1]-1)),-1],1,quantile,probs=0.05)
	c5 <- apply(C[1:n.esp,-1],1,quantile,probs=0.05)
	a95 <- apply(A[1:n.ideg,-1],1,quantile,probs=0.95)
	aa95 <- apply(AA[1:n.odeg,-1],1,quantile,probs=0.95)
	b95 <- apply(B[-(n.dist:(dim(B)[1]-1)),-1],1,quantile,probs=0.95)
	c95 <- apply(C[1:n.esp,-1],1,quantile,probs=0.95)	
	par(mfrow=c(2,2),oma=c(0,0,3,0),mar=c(4,3,1.5,1))
	
	boxplot(as.data.frame(t(A[1:n.ideg,-1])),
	        xaxt="n",
	        xlab="in degree",
	        ylab="proportion of nodes")
	axis(1,seq(1,n.ideg),seq(0,n.ideg-1))
	lines(A[1:n.ideg,1],lwd=2,col=2)
	lines(a5,col="darkgray")
	lines(a95,col="darkgray")
	
	title("Bayesian goodness-of-fit diagnostics",outer=TRUE)
	
	boxplot(as.data.frame(t(AA[1:n.odeg,-1])),
	        xaxt="n",
	        xlab="out degree",
	        ylab="proportion of nodes")
	axis(1,seq(1,n.odeg),seq(0,n.odeg-1))
	lines(AA[1:n.odeg,1],lwd=2,col=2)
	lines(aa5,col="darkgray")
	lines(aa95,col="darkgray")

	boxplot(as.data.frame(t(B[-(n.dist:(dim(B)[1]-1)),-1])),
	        xaxt="n",
	        xlab="minimum geodesic distance",
	        ylab="proportion of dyads")
	axis(1,seq(1,n.dist),labels=c(seq(1,(n.dist-1)),"NR"))
	lines(B[-(n.dist:(dim(B)[1]-1)),1],lwd=2,col=2)
	lines(b5,col="darkgray")
	lines(b95,col="darkgray")
	
	boxplot(as.data.frame(t(C[1:n.esp,-1])),
	        xaxt="n",
	        xlab="edge-wise shared partners",
	        ylab="proportion of edges")
	axis(1,seq(1,n.esp),seq(0,n.esp-1))
	lines(C[1:n.esp,1],lwd=2,col=2)
	lines(c5,col="darkgray")
	lines(c95,col="darkgray")

	out = list(sim.idegree=A[,-1],
	           sim.odegree=AA[,-1],
	           sim.dist=B[,-1],
	           sim.esp=C[,-1],
	           obs.degree=A[,1],
	           obs.dist=B[,1],
	           obs.esp=C[,1],
	           fun=F)
}     
}

