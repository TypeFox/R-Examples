internds <-
function(x){
	#return internal edges and their lengths, labeled by the descendant node
	n <- length(x$tip.label)
	intin<-x$edge[,2]>n
	anc <- x$edge[intin,1]
	time <- as.numeric(branching.times(x))
	names(time) <- (n+1):(n+x$Nnode)
	data.frame(anc=x$edge[intin,1],dec=x$edge[intin,2],len=x$edge.length[intin],time=time[match(anc,names(time))],label=x$edge[intin,2])
	}

