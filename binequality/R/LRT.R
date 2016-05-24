LRT <-
function(dat,fitComb,ID){
	dists<-unique(fitComb[,ID])
	GFs<-rep(NA, nrow(fitComb))
	dfs<-rep(NA, nrow(fitComb))
	ps<-rep(NA, nrow(fitComb))
	for(i in dists){
		use.i<-which(fitComb[,ID]==i)
		fitComb.i<-fitComb[use.i,]
		dat.i<-dat[which(dat[,ID]==i),]
		n.i<-sum(dat.i[,'hb'],na.rm=TRUE)
		nb_ln.i<-dat.i[,'hb']*log(dat.i[,'hb']/n.i)
		sat.loglik.i<-sum(nb_ln.i,na.rm=TRUE)
		B1.i<-nrow(dat.i)-1
		Bplus.i<-length(which(dat.i[,'hb']>0))
		G2.i<- -2*(fitComb.i$logLikelihood-sat.loglik.i)
		df.i<-min(c(B1.i,Bplus.i),na.rm=TRUE)-fitComb.i[,'nparams']
		p.i<-1-pchisq(G2.i,df.i)
		GFs[use.i]<-G2.i
		dfs[use.i]<-df.i
		ps[use.i]<-p.i
	}
	G2<-GFs
	df<-dfs
	p<-ps
	fitComb.out<-data.frame(fitComb,G2,df,p)
	return(fitComb.out)
}
