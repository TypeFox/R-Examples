DetermineWeight_SilClust<-function(List,type=c("data","dist","clusters"),weight=seq(0,1,by=0.01),distmeasure=c("tanimoto","tanimoto"),normalize=FALSE,method=NULL,nrclusters=NULL,names=NULL,nboot=1000,StopRange=FALSE,plottype="new",location=NULL){
	
	type=match.arg(type)
	
	if(is.null(names)){
		names=c()
		for(i in 1:length(List)){
			names=c(names,paste("Method",i,sep=" "))
		}	
	}
	
	CheckDist<-function(Dist,StopRange){
		if(StopRange==FALSE  & !(0<=min(Dist) & max(Dist)<=1)){
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
	else if(type=="dist"){
		OrderNames=rownames(List[[1]])
		for(i in 1:length(List)){
			List[[i]]=List[[i]][OrderNames,OrderNames]
		}
		Dist=List
		Dist=lapply(seq(length(Dist)),function(i) CheckDist(Dist[[i]],StopRange))
		silwidth=lapply(Dist,function(x) cluster::pam(x,nrclusters)$silinfo$widths)
		names(silwidth)=names
	}
	else{
		Dist=lapply(seq(length(List)),function(i) return(List[[i]]$DistM))
		Dist=lapply(seq(length(Dist)),function(i) CheckDist(Dist[[i]],StopRange))
		OrderNames=rownames(Dist[[1]])
		for(i in 1:length(Dist)){
			Dist[[i]]=Dist[[i]][OrderNames,OrderNames]
		}
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
	
	
	namesw=c()
	for(i in 1:length(names)){
		namesw=c(namesw,paste("w_",names[i],sep=""))
	}
	
	labels<-c(namesw,"Observed Statistic","P-Value")
	
	ResultsWeight<-matrix(0,ncol=length(labels),nrow=length(weight))
	colnames(ResultsWeight)=labels
	
	if(is.null(weight)){
		equalweights=1/length(List)
		weight=list(rep(equalweights,length(List)))		
	}
	else if(class(weight)=='list' & length(weight[[1]])!=length(List)){
		stop("Give a weight for each data matrix or specify a sequence of weights")
	}
	else{
		message('The weights are considered to be a sequence, each situation is investigated')
	}
	
	if(class(weight)!="list"){
		condition<-function(l){		
			l=as.numeric(l)
			if( sum(l)==1 ){  #working wit characters since with the numeric values of comb or permutations something goes not the way is should: 0.999999999<0.7+0.3<1??
				#return(row.match(l,t1))
				return(l)
			}
			else(return(0))
		}
		t1=gtools::permutations(n=length(weight),r=length(List),v=as.character(weight),repeats.allowed = TRUE)
		t2=lapply(seq_len(nrow(t1)), function(i) if(sum(as.numeric(t1[i,]))==1) return(as.numeric(t1[i,])) else return(0)) #make this faster: lapply on a list or adapt permutations function itself: first perform combinations under restriction then perform permutations
		t3=sapply(seq(length(t2)),function(i) if(!all(t2[[i]]==0)) return (i) else return(0))
		weight=t2[which(t3!=0)]
	}
	
	if(class(weight)=="list" & "x" %in% weight[[1]]){ #x indicates a free weight
		for(i in 1:length(weight)){
			w=weight[[i]]
			weightsfordata=which(w!="x") #position of the provided weight = position of the data to which the weight is given
			givenweights=as.numeric(w[weightsfordata])
			
			stilltodistribute=1-sum(givenweights)
			
			newweights=seq(stilltodistribute,0,-0.1)
			
			t1=gtools::permutations(n=length(newweights),r=length(List)-length(weightsfordata),v=as.character(newweights),repeats.allowed = TRUE)
			Input1=as.list(seq_len(nrow(t1)))
			Input2=lapply(seq(length(Input1)),function(i) {Input1[[i]][length(Input1[[i]])+1]=stilltodistribute
						return(Input1[[i]])})
			t2=lapply(seq(length(Input2)), FUN=function(i){if(sum(as.numeric(t1[Input2[[i]][1],])+0.00000000000000002775)==Input2[[i]][2]) return(as.numeric(t1[i,])) else return(0)}) #make this faster: lapply on a list or adapt permutations function itself: first perform combinations under restriction then perform permutations
			t3=sapply(seq(length(t2)),function(i) if(!all(t2[[i]]==0)) return (i) else return(0))
			weightsforotherdata=t2[which(t3!=0)]
			
			new=list()
			for(i in 1:length(weightsforotherdata)){
				w1=weightsforotherdata[[i]]
				new[[i]]=rep(0,length(List))
				new[[i]][weightsfordata]=givenweights
				new[[i]][which(new[[i]]==0)]=w1
			}
			
			weight=new
		}
	}
	
	weightedcomb<-function(w,Dist){
		temp=lapply(seq_len(length(Dist)),function(i) w[i]*Dist[[i]])
		temp=Reduce("+",temp)		
	}
	
	DistW=lapply(weight,weightedcomb,Dist)
	
	silwidthW=lapply(DistW,function(x) cluster::pam(x,nrclusters)$silinfo$widths)
	
	
	StatRSq<-function(silwidthW,silwidth,ordernames,names){
		n=length(silwidth)
		
		regressRSq<-function(silwidthweight,silwidth,ordernames,names){
					
			L1=silwidthweight[,3][ordernames]		
			L2W=silwidthweight[,1][ordernames]
			
			regressWW<-stats::lm(L1~L2W)
			RsqWW<-summary(regressWW)$r.squared
			
			regressWX<-function(silw,L1,ordernames){
				L2<-silw[,1][ordernames]
				
				regresWX<-stats::lm(L1~L2)
				RsqWX<-summary(regresWX)$r.squared
				return(RsqWX)
			}
			
			RsqWX<-sapply(silwidth,regressWX,L1=L1,ordernames=ordernames)
			
			Rsq=c(RsqWW,RsqWX)
			return(Rsq)
		}
		
		
		RSqs=lapply(c(1:length(silwidthW)),function(x) regressRSq(silwidthW[[x]],silwidth,ordernames,names))
		
		statfunction<-function(RS){
			stat=0
			xx=0
			xy=0
			for(i in 1:length(RS)){
				if(i==1){
					xx=xx+RS[i]
				}
				else{
					xy=xy+RS[i]
				}
			}
			stat=abs(n*xx-xy)  #check this formula with Nolen
			names(stat)=NULL
			return(stat)
		}
		 
		Stats=sapply(RSqs,statfunction)	
		
	}	
	
	StatRSqObs=StatRSq(silwidthW,silwidth,ordernames=rownames(Dist[[1]]),names)

	#bootstrapping
	statNULL=matrix(0,nrow=length(weight),ncol=nboot)
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
	
		
		DistWNULL=lapply(weight,weightedcomb,DistNULL)
		silwidthWNULL=lapply(DistWNULL,function(x) cluster::pam(x,nrclusters)$silinfo$widths)
		
		statNULL[,i]=StatRSq(silwidthWNULL,silwidthNULL,ordernames=rownames(DistNULL[[1]]),names)
	}
	
	PVals=lapply(c(1:nrow(statNULL)),function(x) (1+sum(abs(statNULL[x,])<=abs(StatRSqObs[x])))/(nboot+1))
	
	ResultsWeight=t(mapply(c,weight,StatRSqObs,PVals))
	colnames(ResultsWeight)=labels
	
	#Choose weight with smallest observed test statistic

	Weight=ResultsWeight[which.min(abs(ResultsWeight[,3]-0)),c(1:length(List))]
	
	plottypein(plottype,location)
	graphics::plot(x=ResultsWeight[,1],y=ResultsWeight[,"Observed Statistic"],xlim=c(0,max(ResultsWeight[,1])),ylim=c(min(ResultsWeight[,"Observed Statistic"]),max(ResultsWeight[,"Observed Statistic"])),xlab="",ylab="Observed Statistic",pch=19,col="black")
	graphics::points(ResultsWeight[which.min(abs(ResultsWeight[,3]-0)),1],ResultsWeight[which.min(abs(ResultsWeight[,3]-0)),"Observed Statistic"],pch=19,col="red")
	graphics::mtext("Weight Combinations", side=1, line=4)
	graphics::axis(1,labels=paste("Optimal weights:", paste(Weight,collapse=", "),sep=" "), at=ResultsWeight[which.min(abs(ResultsWeight[,3]-0)),1],line=2)
	plottypeout(plottype)
	
	plottypein(plottype,location)
	graphics::plot(x=ResultsWeight[,1],y=ResultsWeight[,"P-Value"],xlim=c(0,max(ResultsWeight[,1])),ylim=c(min(ResultsWeight[,"P-Value"]),max(ResultsWeight[,"P-Value"])),xlab="",ylab="P-Value",pch=19,col="black")
	graphics::points(ResultsWeight[which.min(abs(ResultsWeight[,3]-0)),1],ResultsWeight[which.min(abs(ResultsWeight[,3]-0)),"P-Value"],pch=19,col="red")
	graphics::mtext("Weight Combinations", side=1, line=4)
	graphics::axis(1,labels=paste("Optimal weights:", paste(Weight,collapse=", "),sep=" "), at=ResultsWeight[which.min(abs(ResultsWeight[,3]-0)),1],line=2)
	plottypeout(plottype)
	
	
	out=list()
	out[[1]]=ResultsWeight
	out[[2]]=Weight
	names(out)=c("Result","Weight")

	return(out)
	
}






