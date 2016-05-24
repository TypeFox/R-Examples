`simsij` <-
function(nsims=100,n=256, proportion=P2, binsize=1, thrule="ebayesthresh", van=8, fam="DaubLeAsymm", pl=3, 
prior="laplace",vscale="independent", plotstep=FALSE,a=NA,truncate=FALSE,...){


ans <- matrix(0, nrow=nsims, ncol=2)
estmat<-array(0,dim=c(nsims,n,2))
binmat<-matrix(0,nsims,n)

x<-(1:n)/n
y<-proportion(x)

for (i in 1:nsims)	{
	
	ans[i,] <- plotest(l<-hfdenoise(n=n, proportion=proportion, binsize=binsize, thrule=thrule, van=van, fam=fam, pl=pl, 
	prior=prior, vscale=vscale, plotstep=plotstep,truncate=truncate,...))[1:2]
	estmat[i,,1]<-l$fhat
	estmat[i,,2]<-l$fhata
	binmat[i,]<-l$b
	}

return(list(x=x,truep=y,ans=ans,est=estmat,bin=binmat))

}

