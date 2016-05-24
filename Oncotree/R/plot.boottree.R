
"plot.boottree" <-
function(x, minfreq=NULL, minprop=NULL, nboots=NULL, draw.orig=TRUE, draw.consensus=TRUE, 
         fix.nodes=FALSE, ask=(prod(par("mfrow"))<ntrees)&&dev.interactive(),  ...){

  n.non.null <- 3-sum(sapply(list(minfreq, minprop, nboots), is.null)) 
  if (n.non.null > 1) 
  	stop("specify no more than one of 'minfreq', 'minprop', and 'nboots'")
 	if (n.non.null == 0)
 	  minfreq <- 1  #by default draw all trees
  if (!is.null(minfreq))
    nboots <- sum(x$tree.list$Freq >= minfreq) 
  if (!is.null(minprop))
    nboots <- sum(x$tree.list$Freq >= minprop*sum(x$tree.list$Freq))
	nboots <- min(nboots, nrow(x$tree.list)) #in case there are not enough unique trees 
    
	ntrees <- nboots + draw.orig + draw.consensus
  if (ask) {
    op <- par(ask = TRUE)
    on.exit(par(op))
    } 
	
	plotinfo <- build.plot(x$original)
	orig.tree <- list(parent=x$original, level=plotinfo$level, levelgrp=plotinfo$levelgrp, nmut=plotinfo$nmut)

	dots <- list(...)
  cex.sub <- ifelse(is.null(dots$cex),par("cex"),dots$cex) * par("cex.sub") 
  if (fix.nodes){
     if (is.null(dots$node.coords)){
        dots$node.coords <- plot.oncotree(orig.tree, plot=FALSE)
     }
  }
  else {
     if (!is.null(dots$node.coords)){
        dots$node.coords <- NULL
     }
  }
    
  
  if (draw.orig){
#		plot.oncotree(orig.tree,  ...)
		do.call(plot.oncotree, c(x=list(orig.tree), dots))
		mtext("Original Tree", side=1, line=1, cex=cex.sub)
	}
	
	if (draw.consensus){
		child <- x$original$child
		parent.num <- as.numeric(x$consensus)
		parent <- character()
		nmut <- length(child)
		    for (i in 1:nmut){
		        number <- parent.num[i]
		        if (i==1){
		            parent[i] <- ""}
		        else{
		            parent[i] <- child[number]}
		        }
		mostfreq <- list(child=child, parent=parent, parent.num=parent.num)
		plotinfo <- build.plot(mostfreq)
		otree <- list(parent=mostfreq, level=plotinfo$level, levelgrp=plotinfo$levelgrp, nmut=plotinfo$nmut)
		do.call(plot.oncotree, c(x=list(otree), dots))
		mtext("Tree based on most frequent parent", side=1, line=1, cex=cex.sub)
	}
	
	#only unlist the top 'nboots' trees by freq (note: this comes sorted, but we won't use that fact)
	if (nboots > 0){
	  tree.freq <- x$tree.list[order(-x$tree.list$Freq)[1:nboots],]
		treelist <-as.matrix(tree.freq$Tree)
		z <- apply(treelist, 1, function(y) as.vector(as.numeric(unlist(strsplit(y, "\\.")))))
		z2 <- t(z)
		#z2 contains a list of parents for each tree with freq>=minfreq
		nmut <- ncol(z2)
		parent<-character()
		for (p in 1:nrow(z2)){
		    child <- x$original$child
		    parent.num <- z2[p,]
		    parent <- c("",child[parent.num])  #child[0]==character()
		    parent2 <- list(child=child, parent=parent, parent.num=parent.num)
		    keepparent <- which(!is.na(parent2$parent.num))
		    parent4 <- lapply(parent2, function(x)x[keepparent])
		    plotinfo <- build.plot(parent4)
		    otree <- list(parent=parent4, level=plotinfo$level, levelgrp=plotinfo$levelgrp, nmut=plotinfo$nmut)
	    	do.call(plot.oncotree, c(x=list(otree), dots))
#		    if (fix.nodes){
#		       nds <- plot.oncotree(orig.tree, plot=FALSE)
#		       plot.oncotree(otree, node.coords=nds, ...)
#		    }
#		    else {
#		      plot.oncotree(otree,  ...)
#        }
		    freq.count <- tree.freq[p,"Freq"]
		    tt <- paste("Observed Frequency", freq.count, sep=" = " )
		    mtext(tt, side=1, line=1, cex=cex.sub)
    }
  }
}
