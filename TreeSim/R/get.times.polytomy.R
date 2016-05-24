get.times.polytomy<-function(tree){
	nodes<-sort(unique(c(tree$edge)))
	ttype<-(1:length(nodes))*0
	times<-ttype
	ttype[tree$edge[1,1]]<-1
	timespoly<-vector()
	for (j in (tree$edge[1,1]+1):length(nodes)) {
		ttype[j]<-1
		temp<-which(tree$edge[,2]==j)
		ancestor <- tree$edge[temp,1]
		times[j] <- times[ancestor]+tree$edge.length[temp]	
		polytomy <- length(which(tree$edge[,1]==j))
		if (polytomy>2) {temp2<-rep(times[j],polytomy-2)
			timespoly<-c(timespoly,temp2)
		}
	}
	for (j in 1:(tree$edge[1,1]-1)) {
		temp<-which(tree$edge[,2]==j)
		ancestor<-tree$edge[temp,1]
		times[j]<-times[ancestor]+tree$edge.length[temp]	
	}
	times<-c(times,timespoly)
	ttype<-c(ttype,(timespoly*0+1))
	maxt<-max(times)
	times<- -times+maxt
	out<-cbind(times,ttype)
	out
}

