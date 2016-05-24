Hmemory<-function(iicc,training.set){
	nucleotids<-c("A","T","C","G")
	q 	<- iicc$q
	Prob<-iicc$background
	#missing.fun<-iicc$missing.fun
	Entropy<-lapply(seq(1, length(nucleotids), 1), function(i){
					training.set.mes.val<-rbind(as.matrix(training.set),rep(nucleotids[i],ncol(training.set)))
					pmX<-probability(training.set.mes.val, Prob)
					H<-switch(iicc$classentropy, "Shannon"=entropy.Shannon(pmX),"Renyi"=entropy.Renyi(pmX,q))
					H
					})
	
	list(Entropy=Entropy)
	}
