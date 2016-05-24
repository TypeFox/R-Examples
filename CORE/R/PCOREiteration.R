PCOREiteration <-
function(z,maxmark,minscore=0,nprof=1){
chu<-unique(z[,"chrom"])
chc<-rep(0,length(chu))
for(i in 1:length(chu))chc[i]<-sum(z[,"chrom"]==chu[i])
for(ich in 1:length(chu)){
	za<-z[z[,"chrom"]==chu[rev(order(chc))][ich],c("start","end","weight"),drop=F]
	za<-za[order(za[,"start"]),,drop=F]
	winloci<-matrix(ncol=4,nrow=maxmark,
		dimnames=list(NULL,c("chrom","start","end","score")))
	winloci[,"chrom"]<-chu[rev(order(chc))][ich]
	for(i in 1:maxmark){
		if(nrow(za)==0)break
		y<-cbind(c(za[,"end"]+1,za[,"start"]),c(-za[,"weight"],za[,"weight"]))
		y<-y[order(y[,1]),,drop=F]
		cy2<-cumsum(y[,2])
		score<-max(cy2)
		if(score<minscore)break
		stabstart<-y[which.max(cy2),1]
		stabend<-min(y[y[,1]>stabstart&y[,2]<=0,1])-1
		winloci[i,c("start","end","score")]<-c(stabstart,stabend,score)
		if(ich>1)if(((sum(accu[,"score"]>winloci[i,"score"])+i)>=maxmark))break
		za<-za[!(za[,"start"]<=stabstart&za[,"end"]>=stabstart),,drop=F]
	}
	winloci<-winloci[!is.na(winloci[,"score"]),,drop=F]
	if(ich==1)accu<-winloci
	if(ich>1){
		accu<-rbind(accu,winloci)
		accu<-accu[accu[,"score"]>=minscore,,drop=F]
		accu<-accu[order(accu[,"score"],decreasing=T),,
					drop=F][1:min(maxmark,nrow(accu)),,drop=F]
	}
}
accu[,"score"]<-accu[,"score"]/nprof
return(accu)
}
