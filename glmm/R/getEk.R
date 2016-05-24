getEk <-
function(z){
	nrand<-lapply(z,ncol)
	nrandom<-unlist(nrand)
	totnrandom<-sum(nrandom)
	nsigma<-length(z)
	
	Ekmat<-matrix(data=c(0),nrow=nsigma,ncol=totnrandom)
	Ek<-list()

	cumnrand<-c()
	cumnrand[1]<-0
	for(i in 1:length(nrandom)){
		cumnrand[i+1]<-sum(nrandom[1:i])
	}
	
	#each row of Ekmat is a diagonal of one of the Eks

	for(i in 1:nsigma){
		for(j in 1:totnrandom){
			if(cumnrand[i]<j&cumnrand[i+1]>j-1)		Ekmat[i,j]<-1
		}
		Ek[[i]]<-Ekmat[i,]
	}
	Ek
}
