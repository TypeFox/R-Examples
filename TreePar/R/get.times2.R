get.times2<-function(tree){
	nodes<-sort(unique(c(tree$edge)))
	ttype<-(1:length(nodes))*0
	times<-ttype
	ttype[tree$edge[1,1]]<-1
	for (j in (tree$edge[1,1]+1):length(nodes)) {
		ttype[j]<-1
		temp<-which(tree$edge[,2]==j)
		ancestor<-tree$edge[temp,1]
		times[j]<-times[ancestor]+tree$edge.length[temp]	
	}
	for (j in 1:(tree$edge[1,1]-1)) {
		temp<-which(tree$edge[,2]==j)
		ancestor<-tree$edge[temp,1]
		times[j]<-times[ancestor]+tree$edge.length[temp]	
	}
	maxt<-max(times)
	times<- -times+maxt
	out<-cbind(times,ttype)
	out
}
