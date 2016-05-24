
## constructor of DGobj based on a mechanistic model
DGobj.simul.mechanistic=function(sqrtn, size1, size2, theta, beta, M, delta, plots=FALSE){

	## Definition of the spatial grid and the distance matrix
	N=sqrtn^2
	x0=seq(1,sqrt(N),1)
	x=expand.grid(x0,x0)
	D=sqrt(outer(x[,1],x[,1],"-")^2+outer(x[,2],x[,2],"-")^2)

	## Simulation of immigration dates and immigration locations
	nbStrains=length(theta)
	arrivalDates=sample(1:M,nbStrains,replace=TRUE,prob=(((M:1)-1)^2/sum(((M:1)-1)^2)))
	arrivalNcells=rbinom(n=nbStrains,size=N,prob=beta[1]/N)

	## Simulation of the propagation of sub-epidemics
	Y=NULL
	for(i in 1:nbStrains){
		Y[[i]]=matrix(0,N,M)
		Y[[i]][,arrivalDates[i]]=sample(c(rpois(n=arrivalNcells[i],lambda=beta[2]),
			rep(0,N-arrivalNcells[i])),replace=FALSE)
		for(m in (arrivalDates[i]+1):M){
			Y[[i]][,m]=rpois(N,.infection.potential(c(theta[i],delta),D,Y[[i]][,m-1]))
	  	}
	}

	## Final states of all strains at all grid nodes and strain proportions
	Ylast=NULL
	for(i in 1:nbStrains){
		Ylast=cbind(Ylast,Y[[i]][,M])
	}
	pS=t(apply(Ylast,1,function(u) u/sum(u)))

	## Aggregation of the sub-epidemics to obtain the total epidemic
	Yaggreg=0
	for(i in 1:nbStrains){
		Yaggreg=Yaggreg+Y[[i]]
	}

	## Computation of the Growth variables
	Z=log(1+Yaggreg[,M])-log(1+Yaggreg[,M-1])

	## Sampling
	ech=sample(1:N,size=min(size1,sum(Yaggreg[,M]>0)),replace=FALSE,
		prob=0+(Yaggreg[,M]>0))
	frequ=matrix(0,length(ech),nbStrains)
	for(i in 1:length(ech)){
		isols=sample(rep(1:nbStrains,Ylast[ech[i],]),size=min(size2,sum(Ylast[ech[i],])),
			replace=FALSE)
		for(j in 1:nbStrains){
			frequ[i,j]=sum(isols==j)
		}
	}
	effective.sample.size=sum(frequ)

	if(plots){
		## Production of sub-epidemic plots
		par(mfrow=c(nbStrains,M+1),mar=c(0.1,1,1,0.1)*1.5)
		Ymax=0
		for(i in 1:nbStrains){
			Ymax=max(Ymax,max(Y[[i]]))
		}
		for(i in 1:nbStrains){
			for(m in 1:M){
				if(i==1){ main.col=paste("Time",m) } else { main.col="" }
				if(m==1){ main.row=paste("Strain",i) } else { main.row="" }
				plot(x,pch=19,cex=Y[[i]][,m]/Ymax*3,asp=1,xlab="",ylab="",axes=F,
					main=main.col)
				mtext(main.row,side=2)
				box()
			}
			plot(x,cex=pS[,i]*3,pch=19,asp=1,axes=F,xlab="",ylab="",main=paste("Prop. strain ",i))
			box()
		}
	}

	## Construction of the genetic and demographic data sets
	geneticData=cbind(x[ech,],frequ)
	demographicData=cbind(x,Z)
	colnames(geneticData)=c("x1","x2",paste("strain",1:nbStrains,sep=""))
	colnames(demographicData)=c("x1","x2","growth")
	
    return(new("DGobj",demographic=as.matrix(demographicData),genetic=as.matrix(geneticData)))

}

