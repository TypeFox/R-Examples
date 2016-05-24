

## Ranking method
ranking.strains=function(DGobject, bw, nb.mcsimul, plots=FALSE){

    geneticData=DGobject@genetic
    x=DGobject@demographic[,1:2]
    Z=DGobject@demographic[,3]
    
	nbStrains=ncol(geneticData)-2
	N=nrow(x)
	
	## Estimation of strain proportions
	pSEstimes=matrix(0,nrow = N, ncol = nbStrains)
	for (s in 1:nbStrains) {
		for (i in 1:N) {
    		pSEstimes[i,s]=.estimation.piS(geneticData,x[i,1],x[i,2],s,bw)
  		}
	}

	## Estimation of regression coefficients
	regression <- lm(Z ~ -1 + pSEstimes)
	z=regression$coef

	## Assessment of the significance of differences in the coefficients
	zStar=matrix(0,nb.mcsimul,length(z))
	for(i in 1:nb.mcsimul){
		permut=(1:N)[!is.na(pSEstimes[,1])]
		pSEstimesStar=pSEstimes
		pSEstimesStar[permut,]=pSEstimesStar[sample(permut),]
		regressionStar=lm(Z ~ -1 + pSEstimesStar)
		zStar[i,]=regressionStar$coefficients
	}

	## Computation of pairwise t-tests
	pvals=NULL
	pvals.name=NULL
	for(i in 1:(ncol(zStar)-1)){
		for(j in (i+1):ncol(zStar)){
			pvals=c(pvals,mean( zStar[,j]-zStar[,i] >= z[j]-z[i] ))
			pvals.name=c(pvals.name,paste("(",j,"-",i,")>0",sep=""))
			pvals=c(pvals,mean( zStar[,j]-zStar[,i] <= z[j]-z[i] ))
			pvals.name=c(pvals.name,paste("(",j,"-",i,")<0",sep=""))
			pvals=c(pvals,mean( abs(zStar[,j]-zStar[,i]) >= abs(z[j]-z[i]) ))
			pvals.name=c(pvals.name,paste("(",j,"-",i,")neq0",sep=""))
		}
	}
	names(pvals)=pvals.name
		
	if(plots){	
		## Production of plots (ranking plots)
		nrow.plot=ceiling((2+nbStrains)/3)
		par(mfrow=c(nrow.plot,3),mar=c(2,2,1,0.2))
		plot(x,cex=(Z-min(Z))/(max(Z)-min(Z))*3,pch=19,asp=1,axes=F,xlab="",ylab="",
			main="Growth variable")
		plot(x,cex=1,pch=1,asp=1,axes=F,xlab="",ylab="",main="Sampling sites")
		points(geneticData[,1:2],pch=19)
		for(j in 1:nbStrains){
			zseq=seq(min(c(z,zStar)),max(c(z,zStar)),l=20)
			hist(zStar[,j],xlim=range(zseq),breaks=zseq,xlab="",ylab="",main=paste("Coef. strain ",j))
			abline(v=z[j],lwd=2,lty="dashed",col=2)
		}
	}
	
	names(z)=paste("strain",1:nbStrains,sep="")	
	colnames(zStar)=paste("strain",1:nbStrains,sep="")	

	return(list(permutation.estimates=zStar,estimates=z,p.values=pvals))

}


