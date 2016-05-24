surfaceAICMultiPlot<-function(fwd=NULL, bwd=NULL, out=NULL, summ=NULL, cols=NULL, daic=FALSE, ...){

if(!is.null(out)){
	fwd<-bwd<-list()
	for(i in 1:length(out)){
		fwd[[i]]<-out[[i]][[1]]
		bwd[[i]]<-out[[i]][[2]]
		}
	}

aics<-nreg<-list()
if(is.null(summ)){
	summ<-list()
	for(i in 1:length(fwd)){
		summ[[i]]<-surfaceSummary(fwd[[i]], bwd[[i]])		
		}
	}
for(i in 1:length(summ)){
	aics[[i]]<-summ[[i]]$aics
	nreg[[i]]<-summ[[i]]$n_regimes_seq[2,]
	if(daic)aics[[i]]<-aics[[i]]-aics[[i]][1]
}

if(is.null(cols))cols<-rep("black",length(aics))
	plot(NA,type="n",ylim=range(unlist(aics)),xlim=range(unlist(nreg)),...)
	for(i in 1:length(summ)){
		points(nreg[[i]],aics[[i]],type="l",col=cols[i], ...)
		points(nreg[[i]],aics[[i]],bg=cols[i],pch=21, ...)
	}
abline(h=aics[[i]][1],lty="longdash")	
}
