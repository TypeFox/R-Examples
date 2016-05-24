color.tree.plot<-function(out, tree, p.thres=1, evid.thres=0, PorE=1, show.node.label=FALSE, NODE=TRUE, PADJ=NULL,scale=1,col.rank=TRUE,breaks=50,...)
{	
	################################################################################
	#Empty vectors to store width and colors of edges and nodes
	################################################################################
	ecolor <- rep(NA, length(tree$edge.length))
	ewidth <- rep(1, length(tree$edge.length))


#	subout <- out[out$p.val<p.thres & !is.na(out$p.val),]
	subout <- out[!is.na(out$p.val),]

	if(dim(subout)[1]<1)
	{	print("Not adequate subtree evaluations or very low p threshold");
		print("Nothing to be done")
		return()
	}
	edges <- subout$node2
	tree$node.label <- (length(tree$tip.label)+1):(2*length(tree$tip.label)-1)
	

	################################################################################
	# Calc the ratios of expected rates of each subtree based on its best model
	################################################################################

	ratio <- c()

	for(i in 1:length(edges))
	{	r1 <- switch(subout$mod.1r.tot[i],
			subout$Par1.tot[i],
			calc.rate.weib(as.numeric(subout[i,1:2])),
			calc.rate.lnorm(as.numeric(subout[i,1:2])),
			calc.rate.vrat(as.numeric(subout[i,1:2]))
			)
		r2 <- switch(subout$mod.2r.tr2[i],
			subout$Par1.tr2[i],
			calc.rate.weib(as.numeric(subout[i,5:6])),
			calc.rate.lnorm(as.numeric(subout[i,5:6])),
			calc.rate.vrat(as.numeric(subout[i,5:6]))
			)
		ratio <- c(ratio,r2/r1)
	}
	

	################################################################################
	# Control the range of colors to be used
	# This is based on the ratio of relative diversification rates
	################################################################################

	nb <- sum(ratio<1&ratio>0)
	nr <- sum(ratio>1)

	if(any(ratio<0))
	{	print("Warning: Some rate estimates are negative.")
		print("Most likely due to variable rates model with alpha<1")
	}
	
	col <- rep(NA,length(ratio))
	col[ratio<0] <- "black"
	
	if(col.rank)
	{	lb <- 2*nb
		ub <- 3*nb -1
		ord <- rev(order(ratio[ratio<1]))
		col[ratio<1] <- rainbow((4*nb))[lb:ub][ord]

		lr <- 1
		ur <- 2*nr
		ord <- order(ratio[ratio>1])
		col[ratio>1] <- rainbow((5*nr))[lr:ur][ord]
	}
	else
	{	COL <- rainbow(5*breaks)
		red <- rev(COL[1:breaks])
		ini <- floor(2.5*breaks)
		blue <- COL[ini:(ini+breaks-1)]
		val <- log(ratio[ratio>0])
		max <- max(abs(val))
		br <- max/breaks
		pmarg <- exp(seq(0,max,br))
		nmarg <- exp(seq(0,-max,-br))
		
		for(v in 2:length(pmarg))
		{	col[ratio>pmarg[v-1] & ratio<=pmarg[v] & ratio>0] <- red[v-1]
			col[ratio<nmarg[v-1] & ratio>=nmarg[v] & ratio>0] <- blue[v-1]
		}
	}


	################################################################################
	# Identifying nodes/edges to be highlighted based either on p.val or EvidRatio
	################################################################################

	switch(PorE,
	{	# Adjust p-values based on bonferroni or other. No adjustment if PADJ=NULL
		if(is.character(PADJ))
		{	subout$p.val <- p.adjust(subout$p.val, method=PADJ)
		}
		col <- col[subout$p.val<p.thres]
		ratio <- ratio[subout$p.val<p.thres]
		subout <- subout[subout$p.val<p.thres,]
		edges <- subout$node2
		scl <- 1.5*exp(-2*subout$p.val)*scale
	},
	{	col <- col[subout$E>evid.thres]
		ratio <- ratio[subout$E>evid.thres]
		subout <- subout[subout$E>evid.thres,]
		edges <- subout$node2
		scl <- 1.5*exp(-2*subout$p.val)*scale
	}
	)

	################################################################################
	# Plotting the tree with highlighted nodes or edges
	################################################################################

	if(NODE)			# Color on nodes
	{	plot.phylo2(tree, show.node.label=show.node.label,...)
		edgelabels2(edge=match(edges,tree$edge[,2]),pch=19, col=col,cex=2*scl)
	}
	else				# Color on edges
	{	ecolor[match(edges,tree$edge[,2])] <- col
		ecolor[is.na(ecolor)] <- "black"
		ewidth[match(edges,tree$edge[,2])] <- scl
		
		plot.phylo2(tree, edge.color=ecolor, edge.width=ewidth,show.node.label=show.node.label,...)
	}
}
