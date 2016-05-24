add.rootedge <-
function(t1,age){
	t<-t1
	t$edge <- t1$edge
	root<- max(t$edge)+1
	mrca<-t$edge[1,1]
	t$edge <- rbind(c(root,mrca),t$edge)
	edgelength <- age-t$edge.length[1]
	t$edge.length <- c(edgelength,t$edge.length)
	t$tip.label <- paste("t", sample(t$Nnode+1), sep = "")
	t$Nnode<-t$Nnode+1
	t
	}

