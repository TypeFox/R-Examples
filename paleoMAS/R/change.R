change <-
function(x,age,dca=FALSE,meth="euclidean",bin=FALSE,
	roc=FALSE,digits=1)
{
library(vegan)
{
	age<-round(age,digits=digits)
	if(dca==TRUE){
		scores(decorana(x),display="sites")->scores
		diag(as.matrix(vegdist(scores,method=meth))[-1,
			-nrow(x)])->dist.vect
			if(roc==TRUE){
				age[-1]-age[-nrow(x)]->res
				dist.vect/res->r.of.c
				}
			if(roc==FALSE){
				plot(dist.vect,age[-nrow(x)],ylim=c(max(age),
					min(age)),type="l",ylab="Age",
					xlab="DCA distance")
				}
			else{
				plot(r.of.c,age[-nrow(x)],ylim=c(max(age),
					min(age)),type="l",ylab="Age",
					xlab="DCA rate of change")
				}
	}
	else{
		vegdist(x,method=meth,binary=bin)->dist.mat
		diag(as.matrix(dist.mat)[-1,-nrow(x)])->dist.vect
			if(roc==TRUE){
				age[-1]-age[-nrow(x)]->res
				dist.vect/res->r.of.c
				}
			if(roc==FALSE){
				plot(dist.vect,age[-nrow(x)],ylim=c(max(age),
					min(age)),type="l",ylab="Age",
					xlab="Distance")
				}
			else{
				plot(r.of.c,age[-nrow(x)],ylim=c(max(age),
					min(age)),type="l",ylab="Age",
					xlab="Rate of change")
				}
	}
}
	if(roc==TRUE){
	r.of.c<-cbind(age[-length(age)],r.of.c)
	colnames(r.of.c)<-c("Age","Rate of change")
	results<-list(r.of.c,res)
	return(results)
	}
	else{
	dist.vect<-cbind(age[-length(age)],dist.vect)
	colnames(dist.vect)<-c("Age","Distance")	
	return(dist.vect)
	}
}

