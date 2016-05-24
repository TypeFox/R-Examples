plot.blca <-
function(x, which=1L, main="", col1=heat.colors(12), ...){

	logistore<- devAskNewPage()
	devAskNewPage(TRUE)
	
	show<- rep(FALSE, 9)
	show[which]<-TRUE

	if(show[1]){
	 Theta<- x$itemprob
	Tau<- x$classprob
	
	sel<- Tau<1e-6 
	if(any(sel)){
		warning("The following groups have neg`ligible membership probability, and will be omitted from the plotting device: ", "\n", 
		paste(rep("Group", sum(sel)), which(sel), collapse=","), ".")
		Tau<- Tau[!sel]
		Theta<- Theta[!sel,]
	}

	M<- ncol(Theta)
	G<- nrow(Theta)
	xcut<-c(0,cumsum(Tau))
	ycut<-0:M
	
	image.plot(ycut, xcut, t(Theta), axes=FALSE, zlim=c(0,1), ylab="Groups", xlab="Variables", main=main, col=col1)
	#image(ycut, xcut, t(Theta), axes=FALSE, ylab="Groups", xlab="Variables", main=main, ...)
	abline(h=xcut)
	abline(v=ycut)
	
	xlabel<- (xcut[-(G+1)] + xcut[-1])*0.5
	ylabel<- (ycut[-(M+1)] + ycut[-1])*0.5
	axis(1, ylabel, 1:M)
	axis(2, xlabel, 1:G)
	
	mtext("Item Probabilities", 3, 0.25)
	}#Show1 Parameter Mosaicplot
	if(show[2]){ 
		N<- nrow(x$Z)
		G<- ncol(x$Z)
		
		ind1<- N%/%20 + 1
		o1<- order(x$count, decreasing=TRUE)
		o2<- order(x$classprob, decreasing=TRUE)
		
		Z<- matrix( (x$count*x$Z)[o1, o2], N, G )
		
		for(ind2 in 1:ind1){
		
		mmax<- min(ind2*20, N)
		
		ind3<- ((ind2 - 1)*20 + 1):mmax
		
		mosaicplot(Z[ind3, ], col= (1:length(x$classprob)+1)[o2], main=main, las=2, ...)
		mtext("Classification Uncertainty", 3, 0.25)

		}
		}#Show2 Classification Uncertainty
		if(show[3]){
	M<-dim(x$itemprob)[2]
	G<-dim(x$itemprob)[1]
	
	tau<- x$classprob
		 
	if(is.null(dimnames(x$itemprob)[2])){
		caption<-paste(rep("Column",M), 1:M)
		 }
	else caption<- dimnames(x$itemprob)[2][[1]]
	
	for(m in 1:M){
		
		r<- range(x$samples$itemprob[,,m])
		denmaty<-denmatx<- matrix(0,512, G)
		
		for(g in 1:G){
			dentoy<- density(x$samples$itemprob[, g, m], from=r[1], to=r[2])
			denmatx[, g]<- dentoy$x 
			denmaty[, g]<- tau[g]*dentoy$y
		}
				
		plot.mat<- cbind(denmatx[,1], rowSums(denmaty))
		colnames(plot.mat)<- c("Probability", "Density")
		
		plot(plot.mat,  type='l', main=main,...)
		mtext(caption[m], 3, 0.25)
				
		for(g in 1:G) 	lines(denmatx[, g], denmaty[, g], col=g+1, lty=2, lwd=1)
		
		}
	}#Show3 itemprob Gibbs, Boot
		if(show[4]){
		
		G<-dim(x$itemprob)[1]

		denmaty<-denmatx<- matrix(0,512, G)
		
		for(g in 1:G){
			dentoy<- density(x$samples$classprob[,g])
			denmatx[, g]<- dentoy$x 
			denmaty[, g]<- dentoy$y
			}
				
		plot.mat<- cbind(denmatx[,1], denmaty[,1])
		
		colnames(plot.mat)<- c("Probability", "Density")

		plot(plot.mat,  type='n', main=main, xlim=range(denmatx), ylim=c(0, max(denmaty)), ...)
		mtext("Condit. Membership", 3, 0.25)
		for(g in 1:G) 	lines(denmatx[, g], denmaty[, g], col=g+1, lty=2, lwd=1)		
		}#Show4 classprob Gibbs, Boot
			if(show[5]){
	M<-dim(x$parameters$itemprob)[2]
	G<-dim(x$parameters$itemprob)[1]
	
	xseq<-seq(0,1,length.out=1000)
	
		 
	if(is.null(colnames(x$itemprob))){
		caption<-paste(rep("Column",M), 1:M)
		}else caption<- colnames(x$itemprob)
	
	for(m in 1:M){
		
		var1<- x$parameters$itemprob[,m,1]*x$parameters$itemprob[, m, 2]/((x$parameters$itemprob[,m,1]+x$parameters$itemprob[, m, 2])^2*(x$parameters$itemprob[, m, 1]+x$parameters$itemprob[, m, 2] + 1))
		
		r1<- max(min(x$itemprob[,m] - 4*sqrt(var1)),0)
		r2<- min(max(x$itemprob[,m] + 4*sqrt(var1)),1)
		
		xseq<-seq(r1,r2,length.out=1000)
		ymax<-0
		
		mden<- x$classprob[1]*dbeta(xseq, x$parameters$itemprob[1, m, 1], x$parameters$itemprob[1, m, 2] )
			
		for(g in 2:G) mden<- mden + x$classprob[g]*dbeta(xseq, x$parameters$itemprob[g,m,1], x$parameters$itemprob[g,m,2] )
		
		plot.mat<- cbind(xseq, mden)
		
		colnames(plot.mat)<- c("Probability", "Density")
		
		plot(plot.mat, main=main, type='n', ...)
		mtext(caption[m], 3, 0.25)
		for(g in 1:G) 	lines(xseq, x$classprob[g]*dbeta(xseq, x$parameters$itemprob[g, m, 1], x$parameters$itemprob[g, m, 2] ), col=g+1 )
		
		lines(xseq,  mden, lty=2, lwd=0.5)

		}
		} #Show 5 VB ItemProb
	
	if(show[6])
	{
		M<-dim(x$parameters$itemprob)[2]
		G<-dim(x$parameters$itemprob)[1]

		xseq<-seq(0,1,length.out=100)
		ymat<- matrix(0, 100, G)
		
		for(g in 1:G) ymat[,g]<- dbeta(xseq, x$parameters$classprob[g], sum(x$parameters$classprob[-g]))
		
			if(any(ymat==Inf)){
			warning("Density occuring on one point for some membership values")
			ymat[ymat==Inf]<- max(ymat[ymat!=Inf]+1)
			}
		
		plot.mat<- cbind(xseq, ymat[, 1])
		
		colnames(plot.mat)<- c("Probability", "Density")
		
		plot(plot.mat, ylim=range(ymat), type='n', main=main, ...)
		mtext("Conditional Membership", 3, 0.25)
		for(g in 1:G) lines(xseq,  ymat[,g], lty=1, lwd=0.5, col=g+1)
		
		} #Show 6 VB ClassProb


				if(show[7]){
		
		lbplot<- cbind(1:length(x$poststore), x$poststore)
		colnames(lbplot)<- c("Iteration", "Log-Posterior")
		plot(lbplot, type='b')
		mtext("Algorithm Convergence", 3, 0.25)
		}## Show 7 EM diagnostic
			if(show[8]){
		
	M<-dim(x$itemprob)[2]
	G<-dim(x$itemprob)[1]
		 
	if(is.null(dimnames(x$itemprob)[2])){
		caption<-paste(rep("Column",M), 1:M)
		 }
	else caption<- dimnames(x$itemprob)[2][[1]]
	
	for(m in 1:M){
		matplot(x$samples$itemprob[,G:1,m],  type='l', main=main, sub="Parameter Chains", ylab="Item Probability", xlab="Iteration", col=G:1+1,...)
		mtext(caption[m], 3, 0.25)		
		
		}
		
	matplot(x$samples$classprob[,G:1],  type='l', main=main, sub="Parameter Chains", ylab="Class Probability", xlab="Iteration",col= G:1 +1 ,...)
		mtext("Conditional Membership", 3, 0.25)
		} ## Show8, Gibbs Diagnostic
			if(show[9]){
		
		lbplot<- cbind(1:length(x$lbstore), x$lbstore)
		colnames(lbplot)<- c("Iteration", "Lower Bound")
		plot(lbplot, type='b', ...)
		mtext("Algorithm Convergence", 3, 0.25)
		} ##Show 9, VB diagnostic

	devAskNewPage(logistore)
}
