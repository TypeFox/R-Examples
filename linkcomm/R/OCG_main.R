#######################################################
# Functions for generating and plotting OCGs.         #
#                                                     #
# OCG algorithm from Becker et al. (2012).            #
#                                                     #
# Author: Alex T. Kalinka (alex.t.kalinka@gmail.com)  #
#                                                     #
#######################################################


plot.OCG <- function(x, type = "", ...)
	# S3 method for "plot" generic function.
	# x is an "OCG" object.
	{
	switch(type,
		members = plotLinkCommMembers(x, ...),
		graph = plotOCGraph(x, ...)
		)
	}


print.OCG <- function(x, ...)
	# S3 method for generic function "print".
	# x is an "OCG" object.
	{
	cat("   *** Summary ***\n   Number of nodes = ",x$numbers[2],"\n   Number of edges = ",x$numbers[1],"\n   Number of communities = ",x$numbers[3],"\n   Number of nodes in largest cluster = ",x$clustsizes[1],"\n   Modularity = ",x$modularity,"\n   Q = ",x$Q,"\n")
	}


linkcomm2clustnsee <- function(x, file = "temp.cns", network.name = NULL)
	{
	# x is an object of class "linkcomm" or "OCG".
	if(class(x)!="linkcomm" && class(x)!="OCG"){
		stop("input must be a linkcomm or OCG object\n")
		}
	if(is.null(network.name)){
		fcall <- match.call()
		network.name <- as.character(fcall[2])
		}
	num.clusters <- x$numbers[3]
	nums <- data.frame(clust.size = rep(0,num.clusters), multi.class = rep(0,num.clusters))
	for(i in 1:num.clusters){
		nums[i,1] <- x$clustsizes[match(i,names(x$clustsizes))]
		nodes <- getNodesIn(x,clusterids=i)
		numcl <- x$numclusters[match(nodes,names(x$numclusters))]
		nums[i,2] <- length(which(numcl>1))
		}
	cout <- ""
	for(i in 1:nrow(nums)){
		cout <- paste(cout,i,"(",nums[i,1],",",nums[i,2],"), ",sep="")
		}

	# Pre-amble.
	cat("#ClustnSee analysis export\n#Algorithm:",class(x),"\n#Network:",network.name,"\n#Scope:Network\n#Cluster ID (nb nodes in cluster, nb multi-classed nodes in cluster):\n#",cout,"\n\n\n", file=file, sep="")

	# Clusters and their nodes.
	for(i in 1:num.clusters){
		cat(">ClusterID:",i,"||\n", file=file, sep="", append = TRUE)
		nodes <- getNodesIn(x,clusterids=i)
		for(j in 1:length(nodes)){
			cat(nodes[j],"\n", file=file, sep="", append = TRUE)
			}
		cat("\n", file=file, sep="", append = TRUE)
		}

	return(invisible())

	}


read.OCG <- function(file, elfile = NULL, verbose = FALSE, keep.out = FALSE)
	{

	if(verbose){
		cat("Reading OCG data...\n")
		}

	Modularity <- 0
	Q <- 0

	out <- .C("readOCG", as.character(file), Modularity = as.integer(Modularity), Q = as.double(Q))

	if(!is.null(elfile)){
		edgelist <- read.table(file = elfile, header = FALSE)
	}else{
		cat("No edgelist file specified...\n")
		return
		}

	ret <- list()
	nodes <- unique(c(as.character(edgelist[,1]),as.character(edgelist[,2])))
	cs <- read.table("OCG_numclusters.txt", header = FALSE)
	if(ncol(cs) != 1){
		csr <- as.numeric(cs[,2])
		names(csr) <- as.character(cs[,1])
		mono <- setdiff(nodes,names(csr))
		miss <- rep(1,length(mono))
		names(miss) <- mono
		csr <- append(csr,miss)
		csr <- sort(csr, decreasing = TRUE)
	}else{ # All monoclustered.
		csr <- rep(1,length(nodes))
		names(csr) <- nodes
		}		

	ci <- read.table("OCG_clusters.txt", header = FALSE)
	colnames(ci) <- c("node","cluster")

	ig <- graph.edgelist(cbind(as.character(edgelist[,1]),as.character(edgelist[,2])), directed = FALSE) # Creates "igraph" object.

	ret$numbers <- c(nrow(edgelist), length(nodes), max(ci[,2])) # Numbers of edges, nodes, and clusters.
	ret$modularity <- out$Modularity
	ret$Q <- out$Q
	ret$nodeclusters <- ci
	ret$numclusters <- csr
	ret$igraph <- ig
	ret$edgelist <- cbind(as.character(edgelist[,1]),as.character(edgelist[,2]))

	# Extract the node-size of each edge cluster and order largest to smallest.
	ss <- NULL
	unn <- unique(ci[,2])
	for(i in 1:length(unn)){
		if(verbose){
			mes<-paste(c("Extracting cluster sizes... ",floor((i/length(unn))*100),"%"),collapse="")
			cat(mes,"\r")
			flush.console()
			}
		ss[i] <- length(which(ci[,2]==unn[i]))
		}
	names(ss) <- unn
	ss <- sort(ss,decreasing=T)
	ret$clustsizes <- ss

	class(ret) <- "OCG"

	if(!keep.out){
		file.remove("OCG_numclusters.txt")
		file.remove("OCG_clusters.txt")
		}

	if(verbose){
		cat("\n")
		}

	return(ret)

	}


getOCG.clusters <- function(network, init.class.sys = 3, max.class.card = 0, cent.class.sys = 1, min.class = 2, verbose = TRUE, keep.out = FALSE)
	{

	del <- FALSE

	if(class(network) == "data.frame" || class(network) == "matrix"){
		if(ncol(network)==2){
			numnodes = length(unique(c(as.character(network[,1]),as.character(network[,2]))))
			write.table(network,file="OCG_input.txt",quote=FALSE,sep="\t",col.names=FALSE,row.names=FALSE)
			del <- TRUE
		}else{
			cat("input network must have 2 columns\n")
			return(invisible())
			}
	}else if(class(network) != "character"){
		cat("network must be a data frame or matrix with 2 columns, or a valid file name\n")
		return(invisible())
	}else if(file.access(network) == -1){
		stop(cat("\nfile not found: \"",network,"\"\n",sep=""))
	}else{
		nn <- read.table(file = network, header = FALSE, stringsAsFactors = FALSE)
		numnodes <- length(unique(c(as.character(nn[,1]),as.character(nn[,2]))))
		rm(nn)
		}

	if(verbose){
		verb <- 1
	}else{
		verb <- 0
		}

	if(del){
		file <- "OCG_input.txt"
	}else{
		file <- network
		}

	# Variables for C code:
	#FichE - filename of input file.
	#typ <- as.integer(init.class.sys)
	#FuStyl <- as.integer(fusion.method)
	#CardMax <- as.integer(max.class.card)
	#FCS <- as.integer(cent.class.sys)
	#ClCh <- as.ineteger(min.class)
	
	ICS <- init.class.sys
	FM <- 0
	MCC <- max.class.card
	CCS <- cent.class.sys
	MC <- min.class

	success <- 0

	out <- .C("getOCGclusters", as.character(file), as.integer(ICS), as.integer(FM), as.integer(MCC), as.integer(CCS), as.integer(MC), as.integer(numnodes), as.integer(verb), success = as.integer(success))

	if(out$success == 1){
		ocg <- read.OCG(file = "OCG_temp.txt", elfile = file, verbose = verbose, keep.out = keep.out)
		if(!keep.out){
			file.remove("OCG_temp.txt")
			}
		if(del){
			file.remove("OCG_input.txt")
			}
		return(ocg)
	}else{
		if(del){
			file.remove("OCG_input.txt")
			}
		return
		}

	}


which.communities <- function(x, nodes)
	# x is a "linkcomm" or "OCG" object.
	{
	comms <- unique(as.integer(x$nodeclusters[x$nodeclusters[,1]%in%nodes,2]))
	return(comms)
	}



numberEdgesIn <- function(x, clusterids = 1:x$numbers[3], nodes)
	# x is a "linkcomm" or "OCG" object.
	{
	edge.memb <- list()
	# Get node community membership by edges.
	for(i in 1:length(nodes)){
		prog <- paste(c("   Getting node community edge density...",floor((i/length(nodes))*100),"%\r"),collapse="")
		cat(prog)
		flush.console()
		comms <- intersect(which.communities(x,nodes=nodes[i]),clusterids)
		howmany <- NULL
		# How many edges in each community?
		for(j in 1:length(comms)){
			cedg <- x$edgelist[getEdgesIn(x, clusterids = comms[j]),]
			if(is.null(dim(cedg))){
				howmany[j] <- 1
			}else{
				howmany[j] <- length(which(c(as.character(cedg[,1]),as.character(cedg[,2]))==nodes[i]))
				}
			}
		names(howmany) <- comms

		edge.memb[[as.character(nodes[i])]] <- howmany
		}
	cat("\n")

	return(edge.memb)

	}


plotOCGraph <- function(x, clusterids = 1:x$numbers[3], nodes = NULL, pie.local = TRUE, incident = TRUE, layout = layout.fruchterman.reingold, vertex.radius = 0.03, scale.vertices = 0.05, edge.color = "grey", vertex.label.color = "black", vertex.label.cex = 0.8, pal = brewer.pal(7,"Set2"), shownodesin = 0, vlabel = TRUE, random = TRUE, ...)
	# x is an "OCG" object.
	{
	# Make an edgelist based on clusterids or nodes.
	if(is.null(nodes)){
		el <- x$edgelist[getEdgesIn(x, clusterids = clusterids),]
		ig <- graph.edgelist(el,directed = FALSE)
		nodes <- V(ig)$name
	}else{
		if(incident){
			# Only clusters to which named nodes belong.
			clusterids = which.communities(x, nodes = nodes)
			el <- x$edgelist[getEdgesIn(x, clusterids = clusterids),]
			ig <- graph.edgelist(el,directed = FALSE)
			nodes <- V(ig)$name
		}else{
			# Get all clusters to which named nodes and their first-order neighbourhood belong.
			clusterids = which.communities(x, nodes = nodes)
			el <- x$edgelist[getEdgesIn(x, clusterids = clusterids, all = TRUE),]
			nodes <- unique(c(as.character(el[,1]),as.character(el[,2])))
			ig <- graph.edgelist(el,directed = FALSE)
			nodes <- V(ig)$name
			}
		}
	crf <- colorRampPalette(pal,bias=1)
	cols <- crf(length(clusterids))
	if(random){
		cols <- sample(cols,length(clusterids),replace=FALSE)
		}
	names(cols) <- clusterids

	if(shownodesin == 0){
		vnames <- V(ig)$name
	}else{ # Show nodes that belong to more than x number of communities.
		vnames <- V(ig)$name
		inds <- NULL
		for(i in 1:length(vnames)){
			if(x$numclusters[which(names(x$numclusters)==vnames[i])] < shownodesin){
				inds <- append(inds,i)
				}
			}
		vnames[inds] <- ""
		}
	if(vlabel==FALSE){
		vnames = NA
		}

	# Get node community membership by edges.
	if(pie.local){
		edge.memb <- numberEdgesIn(x, clusterids = clusterids, nodes = nodes)
	}else{
		edge.memb <- numberEdgesIn(x, nodes = nodes)
		}

	cat("   Getting node layout...")
	lay <- layout(ig) # vertex x-y positions; will serve as centre points for node pies.
	lay <- layout.norm(lay, xmin=-1, xmax=1, ymin=-1, ymax=1)
	rownames(lay) <- V(ig)$name
	cat("\n")

	node.pies <- .nodePie(edge.memb=edge.memb, layout=lay, nodes=nodes, edges=100, radius=vertex.radius, scale=scale.vertices)
	cat("\n")

	dev.hold(); on.exit(dev.flush())
	# Plot graph.
	plot(ig, layout=lay, vertex.shape="none", vertex.label=NA, vertex.label.dist=1, vertex.label.color=vertex.label.color, edge.color=edge.color, ...)
	labels <- list()
	# Plot node pies and node names.
	for(i in 1:length(node.pies)){
		yp <- NULL
		for(j in 1:length(node.pies[[i]])){
			seg.col <- cols[which(names(cols)==names(edge.memb[[i]])[j])]
			polygon(node.pies[[i]][[j]][,1], node.pies[[i]][[j]][,2], col = seg.col)
			yp <- append(yp, node.pies[[i]][[j]][,2])
			}
		lx <- lay[which(rownames(lay)==names(node.pies[i])),1] + 0.1
		ly <- max(yp) + 0.02 # Highest point of node pie.
		labels[[i]] <- c(lx, ly)
		}
	# Plot node names after nodes so they overlay them.
	for(i in 1:length(labels)){
		text(labels[[i]][1], labels[[i]][2], labels = vnames[which(nodes==names(node.pies[i]))], cex = vertex.label.cex, col = vertex.label.color)
		}
	}



.nodePie <- function(edge.memb, layout, nodes, edges, radius, scale = NULL)
	{
	# Returns x-y coordinates for node pies using their layout coordinates as centre points.
	base_radius <- radius
	node.pies <- list()

	t2xy <- function(t){
		t2p <- 2*pi*t + 90 * pi/180
		list(x = radius * cos(t2p), y = radius * sin(t2p))
		}
	for(i in 1:length(edge.memb)){
		prog <- paste(c("   Constructing node pies...",floor((i/length(nodes))*100),"%\r"),collapse="")
		cat(prog)
		flush.console()

		x <- edge.memb[[i]]
		x <- c(0, cumsum(x)/sum(x))
		dx <- diff(x)
		if(length(dx)==1){
			if(is.nan(dx)){
				x <- c(0,1)
				dx <- 1
				}
			}

		if(!is.null(scale)){
			if(length(dx) > 1){
				radius <- base_radius + base_radius*scale*length(dx)
			}else{
				radius <- base_radius
				}
			}		

		for (j in 1:length(dx)) {
			n <- max(2, floor(edges * dx[j]))
			P <- t2xy(seq.int(x[j], x[j + 1], length.out = n))
			if(n != 100){
				P$x <- c(P$x, 0)
				P$y <- c(P$y, 0)
				}
			P$x <- P$x + layout[i,1]
			P$y <- P$y + layout[i,2]		
			mat <- cbind(P$x,P$y)
			colnames(mat) <- c("x","y")
			node.pies[[as.character(nodes[i])]][[j]] <- mat
			}

		}

	return(node.pies)

	}






