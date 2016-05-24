countGraphs <-
function(nodes,max.indeg) {
	grphs<-1
	for (i in 1:max.indeg) {
		grphs<-grphs+choose(nodes,i)
	}
	return(grphs)
}
