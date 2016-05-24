getPermutation <- function(sam.dat,oo,Np=100,type=1,para){
	cat('\ngetting permutation.')
	rlist<-c()
	Ns<-dim(sam.dat)[1]
	Ns.list<-1:Ns
	clist<-Ns.list[!Ns.list%in%oo$list]
	for(i in clist){
		rlist<-c(rlist,rep(i,sam.dat[i,5]))
	}
	Nsp<-floor(para$res.r*length(clist))
	if(Nsp<=2){return(NA)}
	plist<-c()
	tp<-type
	for(i in 1:Np){
		cat('.')
		tmp<-sample(rlist)[1:Nsp]
		sam.per<-sam.dat[tmp,]
		p1.list<-seq(0.05,0.95,by=0.1)
		SD.min<-99999
		cali.best<-NULL
		id<-rownames(sam.per)[1]
		for(p in p1.list){
			for(model in c(1,2)){
				para$model<-model
				cali<-getSumDist(oo,p,sam.per,type=tp,para=para)
				if(cali$SD<SD.min){
					SD.min<-cali$SD
					cali.best<-cali
					cali.best$type<-tp
					cali.best$sam.p<-p
					cali.best$model<-model
				}
			}
		}
		p2.list<-seq(max(cali.best$sam.p-0.1,0.01),min(cali.best$sam.p+0.1,0.99),by=0.01)
		for(p in p2.list){
			for(model in c(1,2)){
				para$model<-model
				cali<-getSumDist(oo,p,sam.per,type=tp,para=para)
				if(cali$SD<SD.min){
					SD.min<-cali$SD
					cali.best<-cali
					cali.best$type<-tp
					cali.best$sam.p<-p
				}
			}
		}
		plist<-c(plist,cali.best$sam.p)
	}
	cat('\n')
	return(plist)
}