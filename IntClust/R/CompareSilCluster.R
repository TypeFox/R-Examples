CompareSilCluster<-function(List,type=c("data","dist"),distmeasure=c("tanimoto","tanimoto"),normalize=FALSE,method=NULL,nrclusters=NULL,names=NULL,nboot=1000,StopRange=FALSE,plottype="new",location=NULL){
	
	type=match.arg(type)
	
	if(is.null(names)){
		names=c()
		for(i in 1:length(List)){
			names=c(names,paste("Method",i,sep=" "))
		}	
	}
	CheckDist<-function(Dist,StopRange){
		if(StopRange==FALSE & !(0<=min(Dist) & max(Dist)<=1)){
			message("It was detected that a distance matrix had values not between zero and one. Range Normalization was performed to secure this. Put StopRange=TRUE if this was not necessary")
			Dist=Normalization(Dist,method="Range")
		}
		else{
			Dist=Dist
		}
		
	}
	
	
	if(type=="data"){
		OrderNames=rownames(List[[1]])
		for(i in 1:length(List)){
			List[[i]]=List[[i]][OrderNames,]
		}
		Dist=lapply(seq(length(List)),function(i) Distance(List[[i]],distmeasure[i],normalize,method))
		Dist=lapply(seq(length(Dist)),function(i) CheckDist(Dist[[i]],StopRange))
		silwidth=lapply(Dist,function(x) cluster::pam(x,nrclusters)$silinfo$widths)
		names(silwidth)=names
		
	}
	else{
		OrderNames=rownames(List[[1]])
		for(i in 1:length(List)){
			List[[i]]=List[[i]][OrderNames,OrderNames]
		}
		Dist=List
		Dist=lapply(seq(length(Dist)),function(i) CheckDist(Dist[[i]],StopRange))
		silwidth=lapply(Dist,function(x) cluster::pam(x,nrclusters)$silinfo$widths)
		names(silwidth)=names
	}
	
	
	plottypein<-function(plottype,location){
		if(plottype=="pdf" & !(is.null(location))){
			grDevices::pdf(paste(location,".pdf",sep=""))
		}
		if(plottype=="new"){
			grDevices::dev.new()
		}
		if(plottype=="sweave"){
			
		}
	}
	
	plottypeout<-function(plottype){
		if(plottype=="pdf"){
			grDevices::dev.off()
		}
	}
	
	regressioncomb=gtools::permutations(n=length(List),r=2,repeats.allowed=T)
	
	StatRSq<-function(regressioncomb,silwidth,ordernames,names){
		
		regressRSq<-function(x,silwidth,ordernames,names){
			
			i1=x[1]
			i2=x[2]
			L1=silwidth[[i1]][,3][ordernames]
			L2=silwidth[[i2]][,1][ordernames]
			
			regress<-stats::lm(L1~L2)
			Rsq<-summary(regress)$r.squared
			#names(Rsq)=paste("RSquared_",names[i1],names[i2],sep="_")
			return(Rsq)
			
			#paste names on this object!!!
		}
		
		RSqs=apply(regressioncomb,1,function(x) regressRSq(x,silwidth,ordernames,names))
		
		#for (i in 1:nrow(regressioncomb)){
		#	names(RSqs)[i]=paste("RSquared_",names[regressioncomb[i,1]],names[regressioncomb[i,2]],sep="_")
		#}
		
		stat=0
		xx=0
		xy=0
		for(i in 1:nrow(regressioncomb)){
			if(regressioncomb[i,1]==regressioncomb[i,2]){
				xx=xx+RSqs[i]
			}
			else{
				xy=xy+RSqs[i]
			}
		}
		stat=abs(xx-xy)  #check this formula with Nolen
		names(stat)=NULL
		return(stat)
		
	}	
	
	StatRSqObs=StatRSq(regressioncomb,silwidth,ordernames=rownames(Dist[[1]]),names)
	
	
	#bootstrapping
	statNULL=c(1:nboot)
	perm.rowscols <- function (D, n) 
	{
		s <- sample(1:n)
		D=D[s, s]
		return(D)
	}
	
	for(i in 1:nboot){
		set.seed(i)
		DistNULL=Dist
		DistNULL[[1]] <- perm.rowscols(DistNULL[[1]],nrow(DistNULL[[1]]))
		
		silwidthNULL=lapply(DistNULL,function(x) cluster::pam(x,nrclusters)$silinfo$widths)
		
		statNULL[i]=StatRSq(regressioncomb,silwidthNULL,ordernames=rownames(DistNULL[[1]]),names)
		
	}
	
	pval=(sum(abs(statNULL)<=abs(StatRSqObs))+1)/(nboot+1)
	
	plottypein(plottype,location)
	graphics::plot(stats::density(statNULL),type="l",main="The Density of the Statistic under the H0")
	graphics::abline(v=StatRSqObs)
	
	out=list()
	out[[1]]=StatRSqObs
	#out[[2]]=statNULL
	out[[2]]=pval
	names(out)=c("Observed Statistic","P-Value")
	
	return(out)
}
