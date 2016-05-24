DPfitSamples <- function(dd,alpha=0.05,low.thr=0.05,min.peaksize=10,prior,mcmc,nt=FALSE){
	sID<-unique(rownames(dd))
	DD<-c()
	sLabel<-rep('case2',length(sID))
	Pval<-c()
	names(sLabel)<-sID
	paraList<-c()
    para<-NULL
	para$model<-1
	para$is.tri<-FALSE
    foo<-function(id){
        sam.dd <- dd[rownames(dd)==id,]
        tmp<-getDPfit(sam.dd[,8],alpha=alpha,low.thr=low.thr,prior=prior,mcmc=mcmc)
        return(tmp)
    }
    if(nt){
        tmpList=mclapply(sID,foo)
    }else tmpList=lapply(sID,foo)
	for(kk in 1:length(tmpList)){
        tmp<-tmpList[[kk]]
        id<-sID[kk]
        sam.dd<-dd[rownames(dd)==id,]
        if(tmp$model==0){
            sLabel[id] <- 'case0'
            paraList <- c(paraList,list())
        }else if(tmp$model==1){
			sLabel[id]<-'case1'
            PA0<-tmp$PA0
            Pval<-c(Pval,tmp$P)
			paraList<-c(paraList,list(A=tmp$A,mu=tmp$mu,sigma=tmp$sigma))
		}else{
            PA0<-tmp$PA0
            Pval<-c(Pval,tmp$P)
			tt<-table(PA0)
			paraList<-c(paraList,list(c()))
			vv<-names(tt)[which(tt<min.peaksize)]
			for(k in vv)PA0[which(PA0==k)]<-NA
			tt<-sort(unique(as.numeric(na.omit(PA0))))
			for(i in 1:length(tt))PA0[which(PA0==tt[i])]<-i
			sam.dd<-cbind(sam.dd,PA0)
			DD<-rbind(DD,sam.dd)
		}
	}
	names(paraList)<-sID
	return(list(DD=DD,Labels=sLabel,Pval=Pval,par=paraList))
}