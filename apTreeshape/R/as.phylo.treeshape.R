"as.phylo.treeshape" <-
function(x, ...) {

set.height <- function(edge) {

	res=seq(0,0, length=nrow(edge))
	node=seq(0,0,length=nrow(edge))
	current=0
	for (i in 1:nrow(edge)) {
		if (edge[i, 2]<0) {
			res[i]=current-node[-edge[i,1]]+1
			node[-edge[i, 2]]=current+1
			current=current+1
		}
	}
	
	total=1+max(node)
	for (i in 1:nrow(edge)) {
		if (edge[i, 2]>0) {
			res[i]=total-node[-edge[i,1]]
		}
	}

	res

}
	tree=x

	edge<-matrix(NA, nrow=(2*nrow(tree$merge)), ncol=2)
	
	current.node=-2
	current.line=1
	
	tmp1=c(nrow(tree$merge), 1, -1)
	tmp2=c(nrow(tree$merge), 2, -1)
	to.do=c(tmp1, tmp2)

	while(identical(to.do, numeric(0)) == FALSE){
		tmp=list(node=to.do[1],son=to.do[2], ancestor=to.do[3])
		to.do=to.do[-(1:3)]
		
		if (tree$merge[tmp$node, tmp$son] < 0) {
			edge[current.line,]=c(tmp$ancestor, -tree$merge[tmp$node, tmp$son])
			current.line=current.line+1
		} else {
			edge[current.line,]=c(tmp$ancestor, current.node)
			tmp1=c(tree$merge[tmp$node, tmp$son], 1, current.node)
			tmp2=c(tree$merge[tmp$node, tmp$son], 2, current.node)
			to.do= c(tmp1, tmp2, to.do)
			
			current.node=current.node-1
			current.line=current.line+1
		}

	}
	
	edge.length=set.height(edge)
	edge[,1]=as.character(edge[,1])
	edge[,2]=as.character(edge[,2])
	
	res=list(edge=edge, edge.length=edge.length, tip.label=tree$names)
	class(res)='phylo'
	res=old2new.phylo(res)

        res

}

