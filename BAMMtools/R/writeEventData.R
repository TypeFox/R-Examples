# a function takes bammdata object and re-write it back into a tree file and an event csv file

writeEventData <- function(ephy, outtreefile, outeventfile){
  
	if(! "bammdata" %in% class(ephy)) {
		stop("Input has to be a bammdata object.\n");
	}
  
	tree <- as.phylo(ephy);
	write.tree(tree, file = outtreefile);
  
	#get all the nodes in the eventData
	nodes <- c();
	for (i in 1:length(ephy$eventData)) {
		nodes <- c(nodes, ephy$eventData[[i]]$node);
	}
	nodes <- as.integer(unique(nodes));
  
	lr_child <- sapply(nodes, function(x) {
		l <- tree$edge[which(tree$edge[,1] == x), 2];
		if (length(l) == 0) {
			return(c(x, NA));
		} else {
			c(min(ephy$downseq[which(ephy$downseq == l[1]):which(ephy$downseq == ephy$lastvisit[l[1]])]), min(ephy$downseq[which(ephy$downseq == l[2]):which(ephy$downseq == ephy$lastvisit[l[2]])]));
		}
	})
	tree_child <- data.frame(node = nodes, leftchild = tree$tip.label[lr_child[1,]], rightchild = tree$tip.label[lr_child[2,]], stringsAsFactors = FALSE);
  
	eventdata <- data.frame(generation = integer(0), leftchild = character(0), rightchild = character(0), abstime = double(0), lambdainit = double(0), lambdashift = double(0), muinit=double(0), mushift = double(0), stringsAsFactors = FALSE);
	
	for (i in 1:length(ephy$eventData)) {
		t <- ephy$eventData[[i]];
		eventdata <- rbind(eventdata, data.frame(generation = rep(i,dim(t)[1]),
				leftchild = tree_child$leftchild[match(t$node,tree_child$node)],
				rightchild = tree_child$rightchild[match(t$node,tree_child$node)],
				abstime = t$time,
				lambdainit = t$lam1, 
				lambdashift = t$lam2,
				muinit = t$mu1, 
				mushift = t$mu2,
				stringsAsFactors = FALSE
			))
	}
	write.csv(eventdata, file = outeventfile, row.names = FALSE, quote = FALSE);
}