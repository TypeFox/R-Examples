sim.genespeciestree <- function(n,numbsim,lambda,mu,frac=1,age=0){
rho<-frac
#first simulate speciation times
if (age==0) {trees <- sim.bd.taxa(n, numbsim, lambda, mu, frac = rho, complete = FALSE, stochsampling = TRUE)} else {
	trees <- sim.bd.taxa.age(n, numbsim, lambda, mu, frac = rho, age=age,mrca=TRUE)}
statistics <-vector()
for (i in 1:length(trees)){
	brtimes<-vector()
	spectree<-vector()
	coaltree<-vector()
	statisticstemp<-c(1:4)*0
	species<-list()
	specieshelp<-list()
	for (j in 1:n){species<-c(species,list(1))
		specieshelp<-c(specieshelp, list(j))
		}
	x<-c(sort(branching.times(trees[[i]])),10^10)
	for (m in 2:length(x)){
	maxstep<-x[m]-x[m-1]
	merge<-sort(sample(1:length(species),2))
	minspec<-c(min(specieshelp[[merge[1]]]),min(specieshelp[[merge[2]]]))
	spectree<-rbind(spectree,minspec)
	species[merge[1]]<- list(c(species[[merge[1]]],species[[merge[2]]]))
	species<-species[-merge[2]]
	specieshelp[merge[1]]<- list(c(specieshelp[[merge[1]]],specieshelp[[merge[2]]]))
	specieshelp<-specieshelp[-merge[2]]
	for (k in 1:length(species)){
		addtimes<-0
		coaltreetemp<-vector()
		if (length(species[[k]])>1){
			while (length(species[[k]])>1 && addtimes<maxstep){
			numb<-length(species[[k]])
			addtimes<-addtimes+rexp(1,(numb*(numb-1)/2))
			if (addtimes<maxstep){
				brtimes<-c(brtimes,(addtimes+x[m-1]))
				merge<-sort(sample(1:length(species[[k]]),2))
				coaltreetemp<- rbind(coaltreetemp,c(addtimes,min(specieshelp[[k]][merge[1]]),min(specieshelp[[k]][merge[2]])))
				statisticstemp<-stats(species[[k]],merge[1],merge[2],statisticstemp)
				species[[k]][merge[1]]<-species[[k]][merge[1]]+species[[k]][merge[2]]
				species[[k]]<-species[[k]][-merge[2]]
				specieshelp[[k]][merge[1]]<-min(specieshelp[[k]][merge[1]],specieshelp[[k]][merge[2]])
				specieshelp[[k]]<-specieshelp[[k]][-merge[2]]
				}	
			}
		}	
		if (length(coaltreetemp)>0) {
			coaltreetemp<-coaltreetemp[order(coaltreetemp[,1]),]
			coaltree<-rbind(coaltree,coaltreetemp)
		}
	}
	}
	treeeq<-0
	coaltree<-coaltree[,2:3]
	if (sum(spectree==coaltree)==length(spectree)){treeeq<-1}
	gammastat<-gamStat(brtimes,return.list=FALSE)
	gammastatnostruc<-gamStat(sort(branching.times(trees[[i]])),return.list=FALSE)
	statistics <-rbind(statistics,c(statisticstemp,treeeq,gammastat,gammastatnostruc))
}
rownames(statistics)<-c(1:length(statistics[,1]))
colnames(statistics)<-c("Colless","s","Sackin","cherries","matching","gammastruc","gamma")
statistics
}
