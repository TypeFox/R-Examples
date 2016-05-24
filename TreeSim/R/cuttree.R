cuttree<-function(tree,cuttime){
	ages<-age.nodes(tree)
	height<-max(ages)
	finaltreeheight<-height-cuttime
	for (i in length(tree$tip.label):1){
		ages<-age.nodes(tree)
		help<-which(tree$edge[,2]==i)
		parent<-tree$edge[help,1]
		if (ages[i]<cuttime){
			if (ages[parent]<cuttime){
				tree<-drop.tip(tree,i)
			} else { 
				tree$edge.length[help]<-ages[parent]-cuttime
			}
		}
	}
tree
}