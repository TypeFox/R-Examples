get.groups<-function(tree,S,xcut=0){
	if (xcut>0) {groups<-cbind(xcut,S)} else {
	ages<-vector()
	for (i in 1:length(S)){
		index<-which(tree$edge[,2]==i)
		temp<-tree$edge.length[index]
		ages<-c(ages,temp)
	}
	groups<-cbind(ages,S)}
	groups
}