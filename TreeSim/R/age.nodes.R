age.nodes<-function(tree){
	ages<-rep(0,length(tree$edge[,1]))
	for (i in 1:length(tree$edge[,1])){
		new<-tree$edge[i,2]
		old<-tree$edge[i,1]
		ages[new]<-ages[old]+tree$edge.length[i]
	}	
	ages<-max(ages)-ages
	ages
}