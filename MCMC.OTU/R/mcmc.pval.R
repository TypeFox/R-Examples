mcmc.pval <-
function(dat,testlim=0,sided=2,ptype="z"){
	dat=data.frame(dat)
	lim=1/length(dat[,1])
	bps=c()
	for (i in c(1:length(names(dat)))){
		ss=dat[,i]-testlim
		mm=mean(ss)
		if(ptype=="z"){
			zs=mean(ss)/sd(ss)
			tst=sided*(1-pnorm(abs(zs)))	
		}
		else {
			if (sided==2) tst=sum(mm*ss<0) 
			if (sided==1) tst=sum(ss<0)
			if (tst==0) tst=sided*lim else tst=sided*tst/length(ss)
#		tst=as.numeric(tst)
		}
		bps=append(bps,tst)
	}
	return(bps)
}
