DetMCDLTS<-function(x,y,intercept=1,alpha=0.75,h=NULL,scale_est="scaleTau2"){
	S2<-DetLTS_raw_dmcd(x=x,y=y,h=h,alpha=alpha,scale_est=scale_est)
	if(length(S2$SubsetSize)==1){
		S1<-ltscheckout(x=x,y=y,inbest=S2$Subset[[1]],h=S2$SubsetSize[1],intercept=intercept,alpha=alpha,use.correction=FALSE,objfct=S2$Objective[1])
		return(S1)
	}
	S1<-vector("list",length(S2$Subset))
	names(S1)<-names(S2$Subset)
	for(i in 1:length(which(S2$SubsetSize<nrow(x))))	S1[[i]]<-ltscheckout(x=x,y=y,inbest=S2$Subset[[i]],h=S2$SubsetSize[i],intercept=intercept,alpha=alpha,use.correction=FALSE,objfct=S2$Objective[i])
	if(max(S2$SubsetSize)==nrow(x))	S1[[i+1]]<-ordreg(x=x,y=y,intercept=intercept)
	S1
}
