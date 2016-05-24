##########################################################
# Function(s) to extract link communitites from		 #
#     directed, undirected, weighted, unweighted, or     #
#     bipartite networks.				 #
# The main algorithms are C++ functions available in     #
#     the R package source.				 #
# 							 #
# Author: Alex T. Kalinka (alex.t.kalinka@gmail.com)	 #
#                                                        #
# See:							 #
#							 #
# Ahn et al. (2010). Link communities reveal multiscale  #
#   complexity in networks. Nature 466:761-765.		 #
#							 #
# Kalinka & Tomancak (2011) linkcomm: an R package for   #
#   the generation, visualization, and analysis of link  # 
#   communities in networks of arbitrary size and type.  #
#   Bioinformatics 27:2011-2012.			 #
#							 #
##########################################################


.onAttach <- function(lib, pkg) 
	{
	packageStartupMessage("\nWelcome to linkcomm version ",packageDescription(pkg)$Version,"\n\nFor a step-by-step guide to using linkcomm functions:\n   > vignette(topic = \"linkcomm\", package = \"linkcomm\")\nTo run an interactive demo:\n   > demo(topic = \"linkcomm\", package = \"linkcomm\")\nTo cite, see:\n   > citation(\"linkcomm\")\nNOTE: To use linkcomm, you require read and write permissions in the current directory (see: help(\"getwd\"), help(\"setwd\"))\n")
	}


integer.edgelist <- function(network)
	# Returns an edge list with integer numbers replacing elements of the network.
	{
	if(!is.character(network)){
		cn <- cbind(as.character(network[,1]),as.character(network[,2]))
	}else{
		cn <- network
		}
	nodes <- unique(as.character(t(cn)))
	ids <- seq(nodes)
	names(ids) <- nodes
	g <- matrix(ids[t(cn)],nrow(cn),ncol(cn),byrow=TRUE)
	ret <- list()
	ret$edges <- g
	ret$nodes <- ids
	return(ret)
	}


edge.duplicates <- function(network, verbose = TRUE)
	# Finds and removes loops, duplicate edges, and bi-directional edges.
	{
	xx <- cbind(as.character(network[,1]),as.character(network[,2]))
	edges <- integer.edgelist(network)$edges
	ne <- nrow(edges)
	loops <- rep(0,ne)
	dups <- rep(0,ne)
	
	out <- .C("edgeDuplicates",as.integer(edges[,1]),as.integer(edges[,2]),as.integer(ne), loops = as.integer(loops), dups = as.integer(dups), as.logical(verbose))

	if(verbose){cat("\n")}

	loops <- which(out$loops == 1)
	dups <- which(out$dups == 1)
	inds <- unique(c(loops,dups))
	ret <- list()
	ret$inds <- inds
	if(length(inds)>0){
		ret$edges <- xx[-inds,]
	}else{
		ret$edges <- xx
		}
	if(verbose){
		if(length(loops)>0){
			cat("   Found and removed ",length(loops)," loop(s)\n",sep="")
			}
		if(length(dups)>0){
			cat("   Found and removed ",length(dups)," duplicate edge(s)\n",sep="")
			}
		}
	return(ret)
	}


getLinkCommunities <- function(network, hcmethod = "average", use.all.edges = FALSE, edglim = 10^4, directed = FALSE, dirweight = 0.5, bipartite = FALSE, dist = NULL, plot = TRUE, check.duplicates = TRUE, removetrivial = TRUE, verbose = TRUE) 
	# network is an edge list. Nodes can be ASCII names or integers, but are always treated as character names in R.
	# If plot is true (default), a dendrogram and partition density score as a function of dendrogram height are plotted side-by-side.
	# When there are more than "edglim" edges, hierarchical clustering is carried out via temporary files written to disk using compiled C++ code.
	{

	if(is.character(network) && !is.matrix(network)){
		if(file.access(network) == -1){
			stop(cat("\nfile not found: \"",network,"\"\n",sep=""))
		}else{
			network <- read.table(file = network, header = FALSE)
			}
		}
	x <- network
	rm(network)

	if(ncol(x)==3){
		wt <- as.numeric(as.character(x[,3]))
		if(length(which(is.na(wt)==TRUE))>0){
			stop("\nedge weights must be numerical values\n")
			}
		x <- cbind(as.character(x[,1]),as.character(x[,2]))
	}else if(ncol(x)==2){
		x <- cbind(as.character(x[,1]),as.character(x[,2]))
		wt <- NULL
	}else{
		stop("\ninput data must be an edge list with 2 or 3 columns\n")
		}

	if(check.duplicates){
		dret <- edge.duplicates(x, verbose = verbose)
		x <- dret$edges
		if(!is.null(wt)){
			if(length(dret$inds) > 0){
				wt <- wt[-dret$inds]
				}
			}
		rm(dret)
		}

	el <- x # Modified edge list returned to user.
	rm(x)
	len <- nrow(el) # Number of edges.
	nnodes <- length(unique(c(as.character(el[,1]),as.character(el[,2])))) # Number of nodes.

	intel <- integer.edgelist(el) # Edges with numerical node IDs.
	edges <- intel$edges
	node.names <- names(intel$nodes)
	numnodes <- length(node.names)

	if(bipartite){
		# Check that network is bipartite.
		big <- graph.edgelist(as.matrix(el), directed = directed)
		bip.test <- bipartite.mapping(big)
		if(!bip.test$res){
			stop("\nnetwork is not bi-partite; change bipartite argument to FALSE\n")
			}
		bip <- rep(1,length(bip.test$type))
		bip[which(bip.test$type==FALSE)] <- 0
		names(bip) <- V(big)$name
		bip <- bip[match(node.names, names(bip))]
		rm(big, bip.test)
	}else{
		bip <- 0
		}

	rm(intel)

	# Switch depending on size of network.
	if(len <= edglim){
		disk <- FALSE
		if(is.null(dist)){
			emptyvec <- rep(1,(len*(len-1))/2)
			if(!is.null(wt)){ weighted <- TRUE}else{ wt <- 0; weighted <- FALSE}
			if(!use.all.edges){
				dissvec <- .C("getEdgeSimilarities",as.integer(edges[,1]),as.integer(edges[,2]),as.integer(len),rowlen=integer(1),weights=as.double(wt),as.logical(directed),as.double(dirweight),as.logical(weighted),as.logical(disk), dissvec = as.double(emptyvec), as.logical(bipartite), as.logical(verbose))$dissvec
			}else{
				dissvec <- .C("getEdgeSimilarities_all",as.integer(edges[,1]),as.integer(edges[,2]),as.integer(len),as.integer(numnodes),rowlen=integer(1),weights=as.double(wt),as.logical(FALSE),as.double(dirweight),as.logical(weighted),as.logical(disk), dissvec = as.double(emptyvec), as.logical(bipartite), as.logical(verbose))$dissvec
				}
			distmatrix <- matrix(1,len,len)
			distmatrix[lower.tri(distmatrix)] <- dissvec
			colnames(distmatrix) <- 1:len
			rownames(distmatrix) <- 1:len
			distobj <- as.dist(distmatrix) # Convert into 'dist' object for hclust.
			rm(distmatrix)
		}else{	
			# Did the user provide an adequate distance matrix?
			if(class(dist) != "dist"){
				stop("\ndistance matrix must be of class \"dist\" (see ?as.dist)\n")
			}else if(attr(dist,which="Size") != len){
				stop("\ndistance matrix size must equal the number of edges in the input network\n")
			}else if(length(dist) != (len*(len-1))/2){
				stop("\ndistance matrix must be the lower triangular matrix of a square matrix\n")
				}
			distobj <- dist
			}
		if(verbose){
			cat("\n   Hierarchical clustering of edges...")
			}
		#if(hcmethod=="energy"){
		#	hcedges <- energy.hclust(distobj)
		#}else{
		#	hcedges <- hclust(distobj, method = hcmethod)
		#	}
		hcedges <- hclust(distobj, method = hcmethod)
		hcedges$order <- rev(hcedges$order)
		rm(distobj)
		if(verbose){cat("\n")}
	}else{
		disk <- TRUE
		if(!is.null(wt)){ weighted <- TRUE}else{ wt <- 0; weighted <- FALSE}
		if(!use.all.edges){
			rowlen <- .C("getEdgeSimilarities",as.integer(edges[,1]),as.integer(edges[,2]),as.integer(len),rowlen=integer(len-1),weights=as.double(wt),as.logical(directed),as.double(dirweight),as.logical(weighted),as.logical(disk), dissvec = double(1), as.logical(bipartite), as.logical(verbose))$rowlen
		}else{
			rowlen <- .C("getEdgeSimilarities_all",as.integer(edges[,1]),as.integer(edges[,2]),as.integer(len),as.integer(numnodes),rowlen=integer(len-1),weights=as.double(wt),as.logical(FALSE),as.double(dirweight),as.logical(weighted),as.logical(disk), dissvec = double(1), as.logical(bipartite), as.logical(verbose))$rowlen
			}
		if(verbose){cat("\n")}
		hcobj <- .C("hclustLinkComm",as.integer(len),as.integer(rowlen),heights = single(len-1),hca = integer(len-1),hcb = integer(len-1), as.logical(verbose))
		if(verbose){cat("\n")}
		hcedges<-list()
		hcedges$merge <- cbind(hcobj$hca, hcobj$hcb)
		hcedges$height <- hcobj$heights

		hcedges$order <- .C("hclustPlotOrder",as.integer(len),as.integer(hcobj$hca),as.integer(hcobj$hcb),order=integer(len))$order
		hcedges$order <- rev(hcedges$order)
		hcedges$method <- "single"
		class(hcedges) <- "hclust"

		}

	# Calculate link densities, cut the tree, and extract optimal clusters.

	hh <- unique(round(hcedges$height, digits = 5)) # Round to 5 digits to prevent numerical instability affecting community formation.
	countClusters <- function(x,ht){return(length(which(ht==x)))}
	clusnums <- sapply(hh, countClusters, ht = round(hcedges$height, digits = 5)) # Number of clusters at each height.
	numcl <- length(clusnums)

	ldlist <- .C("getLinkDensities",as.integer(hcedges$merge[,1]), as.integer(hcedges$merge[,2]), as.integer(edges[,1]), as.integer(edges[,2]), as.integer(len), as.integer(clusnums), as.integer(numcl), pdens = double(length(hh)), heights = as.double(hh), pdmax = double(1), csize = integer(1), as.logical(removetrivial), as.logical(bipartite), as.integer(bip), as.logical(verbose))

	pdens <- c(0,ldlist$pdens)
	heights <- c(0,hh)
	pdmax <- ldlist$pdmax
	csize <- ldlist$csize

	if(csize == 0){
		stop("\nno clusters were found in this network; maybe try a larger network\n")
		}

	if(verbose){
		cat("\n   Maximum partition density = ",max(pdens),"\n")
		}

	# Read in optimal clusters from a file.
	clus <- list()
	for(i in 1:csize){
		if(verbose){
			mes<-paste(c("   Finishing up...1/4... ",floor((i/csize)*100),"%"),collapse="")
			cat(mes,"\r")
			flush.console()
			}
		clus[[i]] <- scan(file = "linkcomm_clusters.txt", nlines = 1, skip = i-1, quiet = TRUE)
		}

	file.remove("linkcomm_clusters.txt")

	# Extract nodes for each edge cluster.
	ecn <- data.frame()
	ee <- data.frame()
	lclus <- length(clus)
	for(i in 1:lclus){
		if(verbose){
			mes<-paste(c("   Finishing up...2/4... ",floor((i/lclus)*100),"%"),collapse="")
			cat(mes,"\r")
			flush.console()
			}
		ee <- rbind(ee,cbind(el[clus[[i]],],i))
		nodes <- node.names[unique(c(edges[clus[[i]],]))]
		both <- cbind(nodes,rep(i,length(nodes)))
		ecn <- rbind(ecn,both)
		}
	colnames(ecn) <- c("node","cluster")
	colnames(ee) <- c("node1","node2","cluster")

	# Extract the node-size of each edge cluster and order largest to smallest.
	ss <- NULL
	unn <- unique(ecn[,2])
	lun <- length(unn)
	for(i in 1:length(unn)){
		if(verbose){
			mes<-paste(c("   Finishing up...3/4... ",floor((i/lun)*100),"%"),collapse="")
			cat(mes,"\r")
			flush.console()
			}
		ss[i] <- length(which(ecn[,2]==unn[i]))
		}
	names(ss) <- unn
	ss <- sort(ss,decreasing=T)

	# Extract the number of edge clusters that each node belongs to.
	unn <- unique(ecn[,1])

	iecn <- as.integer(as.factor(ecn[,1]))
	iunn <- unique(iecn)
	lunn <- length(iunn)
	nrows <- nrow(ecn)

	oo <- rep(0,lunn)

	oo <- .C("getNumClusters", as.integer(iunn), as.integer(iecn), counts = as.integer(oo), as.integer(lunn), as.integer(nrows), as.logical(verbose))$counts

	names(oo) <- unn

	if(verbose){cat("\n")}

	pdplot <- cbind(heights,pdens)

	# Add nodeclusters of size 0.
	missnames <- setdiff(node.names,names(oo))
	m <- rep(0,length(missnames))
	names(m) <- missnames
	oo <- append(oo,m)

	all <- list()

	all$numbers <- c(len,nnodes,length(clus)) # Number of edges, nodes, and clusters.
	all$hclust <- hcedges # Return the 'hclust' object. To plot the dendrogram: 'plot(lcobj$hclust,hang=-1)'
	all$pdmax <- pdmax # Partition density maximum height.
	all$pdens <- pdplot # Add data for plotting Partition Density as a function of dendrogram height.
	all$nodeclusters <- ecn # n*2 character matrix of node names and the cluster ID they belong to.
	all$clusters <- clus # Clusters of edge IDs arranged as a list of lists.
	all$edges <- ee # Edges and the clusters they belong to, arranged so we can easily put them into an edge attribute file for Cytoscape.
	all$numclusters <- sort(oo,decreasing=TRUE) # The number of clusters that each node belongs to (named vector where the names are node names).
	all$clustsizes <- ss # Cluster sizes sorted largest to smallest (named vector where names are cluster IDs).
	all$igraph <- graph.edgelist(el, directed = directed) # igraph graph.
	all$edgelist <- el # Edge list.
	all$directed <- directed # Logical indicating if graph is directed or not.
	all$bipartite <- bipartite # Logical indicating if graph is bipartite or not.

	class(all) <- "linkcomm"

	if(plot){
		if(verbose){
			cat("   Plotting...\n")
			}
		if(len < 1500){ # Will be slow to plot dendrograms for large networks.
			if(len < 500){
				all <- plot(all, type="summary", droptrivial = removetrivial, verbose = verbose)
			}else{ # Slow to reverse order of large dendrograms.
				all <- plot(all, type="summary", right = FALSE, droptrivial = removetrivial, verbose = verbose)
				}
		}else if(len <= edglim){
			oldpar <- par(no.readonly = TRUE)
			par(mfrow=c(1,2), mar=c(5.1,4.1,4.1,2.1))
			plot(hcedges,hang=-1,labels=FALSE)
			abline(pdmax,0,col='red',lty=2)
			plot(pdens,heights,type='n',xlab='Partition Density',ylab='Height')
			lines(pdens,heights,col='blue',lwd=2)
			abline(pdmax,0,col='red',lty=2)
			par(oldpar)
		}else{
			plot(heights,pdens,type='n',xlab='Height',ylab='Partition Density')
			lines(heights,pdens,col='blue',lwd=2)
			abline(v = pdmax,col='red',lwd=2)
			}
		}

	return(all)
	
	}



