#####################################################################################################
# Miscellaneous functions for linkcomm.						                    #
#                                                                                                   #
# Author: Alex T. Kalinka (alex.t.kalinka@gmail.com)                                                #
#                                                                                                   #
# See: 												    #
#												    #
# Ahn et al. (2010). Link communities reveal multiscale complexity in networks. Nature 466:761-765. #
#												    #
# Kalinka & Tomancak (2011) linkcomm: an R package for  					    #
#   the generation, visualization, and analysis of link 					    #
#   communities in networks of arbitrary size and type. 					    #
#   Bioinformatics 27:2011-2012.								    #
#                                                                                                   #
#####################################################################################################


print.linkcomm <- function(x, ...)
	# S3 method for generic function "print".
	# x is a "linkcomm" object.
	{
	cat("   *** Summary ***\n   Number of nodes = ",x$numbers[2],"\n   Number of edges = ",x$numbers[1],"\n   Number of communities = ",x$numbers[3],"\n   Maximum partition density = ",max(x$pdens[,2]),"\n   Number of nodes in largest cluster = ",x$clustsizes[1],"\n   Directed: ",x$directed, "\n   Bi-partite: ",x$bipartite, "\n   Hclust method: ",x$hclust$method,"\n")
	}


linkcomm2cytoscape <- function(x, interaction = "pp", ea= "temp.ea")
	# Writes a Cytoscape edge attribute file to the current directory so that edges are assigned to their respective communities by cluster ID.
	# x is a "linkcomm" object.
	{
	app <- FALSE
	eafunc<-function(y){cat(paste(as.character(y[1])," (",as.character(interaction),") ",as.character(y[2])," = ",y[3],"\n",sep=""),file=ea,append=app); app <<- TRUE}
	out<-apply(x$edges,1,eafunc)
	}

	
getCommunityMatrix <- function(x, nodes = head(names(x$numclusters),20))
	# Returns a binary matrix with nodes as rows, communities as columns, and unit entries indicating membership in a community.
	# x is a "linkcomm" object.
	{
	comms <- unique(x$nodeclusters[as.character(x$nodeclusters[,1])%in%nodes,2]) # Community (cluster) IDs.
	commatrix <- matrix(0,length(nodes),length(comms))
	for(i in 1:length(nodes)){
		cl <- unique(x$nodeclusters[as.character(x$nodeclusters[,1])%in%nodes[i],2])
		minds <- match(cl,comms)
		commatrix[i,minds] <- 1
		}
	return(commatrix)
	}


getNestedHierarchies <- function(x, clusid = 1, verbose = TRUE, plot = TRUE, ids = FALSE)
	# Returns clusters that a user-defined cluster is entirely nested in.
	# x is a "linkcomm" object.
	{
	nestedin <- list()
	comms <- unique(x$nodeclusters[,2])
	clus <- x$nodeclusters[which(x$nodeclusters[,2]==clusid),1]
	for(i in 1:length(comms)){
		if(i==clusid){next}
		assign("name",as.character(i))
		hc <- x$nodeclusters[which(x$nodeclusters[,2]==comms[i]),1]
		if(length(which(clus%in%hc==TRUE))==length(clus)){
			nestedin[[name]] <- hc
			nestid <- i
			}
		}
	if(length(nestedin)==0){
		out <- paste("Cluster ",clusid," is not entirely nested in any other clusters.\n",sep="")
		if(verbose){cat(out)}
	}else{
		if(plot){ # Plots the first nested cluster.
			plotLinkCommGraph(x, clusterids = c(clusid,nestid))
			}
		if(!ids){
			return(nestedin)
		}else{
			return(nestid)
			}
		}
	}


getAllNestedComm <- function(x, verbose = FALSE, plot = FALSE)
	# Returns a list of cluster IDs of clusters that each cluster is entirely nested within.
	# x is a "linkcomm" object.
	{
	nests <- list()
	for(i in 1:length(x$clusters[1:x$numbers[3]])){
		out <- paste(c("   Completed... ",floor((i/length(x$clusters[1:x$numbers[3]]))*100),"%"),collapse="")
		cat(out,"\r")
		flush.console()
		temp <- getNestedHierarchies(x, clusid = i, verbose = verbose, plot = plot, ids = TRUE)
		if(length(temp)!=0){
			assign("name",as.character(i))
			nests[[name]] <- temp
			}
		}
	cat("\n")
	
	return(nests)

	}


getCommunityCentrality <- function(x, nodes = names(x$numclusters), type = "commweight", normalise = TRUE)
	# Returns a measure of community centrality for nodes in the network.
	# x is a "linkcomm" object.
	{
	clusterids <- unique(x$nodeclusters[x$nodeclusters[,1]%in%nodes,2])
	if(type == "commconn"){
		cc <- getCommunityConnectedness(x, clusterids = clusterids, normalise = normalise)
		cd <- rep(0,length(nodes))
		for(i in 1:length(nodes)){
			out <- paste(c("   Completed... ",floor((i/length(nodes))*100),"%"),collapse="")
			cat(out,"\r")
			flush.console()
			for(j in 1:length(cc)){
				ee <- x$edgelist[x$clusters[[clusterids[j]]],]
				nn <- length(unique(c(which(ee[,1]==nodes[i]),which(ee[,2]==nodes[i]))))
				cd[i] <- sum(cd[i], nn*cc[j])
				}
			}
		names(cd) <- nodes
		cat("\n")
		return(cd)

	}else if(type == "commweight"){
		cd <- NULL
		for(i in 1:length(nodes)){
			out <- paste(c("   Completed... ",floor((i/length(nodes))*100),"%"),collapse="")
			cat(out,"\r")
			flush.console()

			cids <- as.integer(unique(x$nodeclusters[x$nodeclusters[,1]%in%nodes[i],2]))

			if(length(cids) > 1){
				cR <- 1-getClusterRelatedness(x, clusterids = cids, cluster = FALSE, verbose = FALSE)
				summed <- 0
				for(j in 1:length(cids)){
					inds <- .getUpperTriIndices(length(cids), which = j)
					summed <- summed + (1-mean(cR[inds]))
					}
				cd[i] <- 1 + summed
			}else if(length(cids)==1){
				cd[i] <- 1
			}else{
				cd[i] <- 0
				}			

			}
		names(cd) <- nodes
		cat("\n")
		return(cd)
		
		}

	}


.getUpperTriIndices <- function(numedg, which = 1)
	{
	# Returns indices for an edge in the upper triangular matrix represented as a vector.
	rows <- NULL
	cols <- NULL
	k <- numedg-1
	for(i in 1:(numedg-1)){
		rows <- append(rows, rep(i,k))
		k <- k-1
		cols <- append(cols, (i+1):numedg)
		}
	ret <- union(which(rows==which),which(cols==which))
	return(ret)
	}


corLinkcommCentrality <- function(x, centrality = "degree", type = "commweight", method = "spearman", plot = TRUE, pch = 20, ...)
	# Returns correlation coefficient between the number of communitites that each node belongs to, and a user-chosen measure of node centrality.
	# x is a "linkcomm" object.
	{
	
	if(centrality == "degree"){
		cn <- degree(x$igraph, v=V(x$igraph))
	}else if(centrality == "betweenness"){
		cn <- betweenness(x$igraph, v=V(x$igraph))
	}else if(centrality == "closeness"){
		cn <- closeness(x$igraph, vids=V(x$igraph))
	}else if(centrality == "constraint"){
		cn <- constraint(x$igraph,nodes=V(x$igraph))
	}

	cc <- getCommunityCentrality(x, type = type)
	names(cn) <- V(x$igraph)$name

	cc <- cc[sort(names(cc))]
	cn <- cn[sort(names(cn))]

	if(plot){
		fit <- summary(lm(cn~cc))
		plot(cc,cn,xlab="Community Centrality",ylab = centrality, pch = pch, ...)
		abline(fit$coefficients[1],fit$coefficients[2])
		}

	return(cor(cc,cn,method=method))

	}


orderCommunities <- function(x, clusterids = 1:x$numbers[3], verbose = TRUE)
	# Puts communities in dendrogram order (useful for plotting graphs).
	# x is a "linkcomm" object.
	{
	clusters <- x$clusters[clusterids]
	dend <- x$hclust$order
	miss <- setdiff(dend,unlist(clusters))
	if(length(miss) > 0){
		minds <- match(miss,dend)
		dend <- dend[-minds]
		}
	id <- dend[1]
	ordered <- list()
	clusids <- NULL
	for(i in 1:length(clusters)){
		if(verbose){
			mes <- paste(c("   Ordering communities according to dendrogram...",floor(i/(length(clusters))*100),"%"),collapse="")
			cat(mes,'\r')
			flush.console()
			}
		for(j in 1:length(clusters)){
			if(!is.na(match(id, clusters[[j]]))){
				ordered[[i]] <- clusters[[j]]
				dend <- dend[-match(clusters[[j]], dend)]
				id <- dend[1]
				clusids[i] <- clusterids[j]
				break
				}
			}
		}

	if(verbose){cat("\n")}

	new<-list()
	new$ordered <- ordered
	new$clusids <- clusids

	return(new)

	}


LinkDensities <- function(x, clusterids = 1:x$numbers[3])
	# Returns the link densities of the communities in "x".
	# x is a "linkcomm" object.
	{
	clusters <- x$clusters[clusterids]
	# Calculate LD for each community.
	ld = NULL
	for(i in 1:length(clusters)){
		edges = length(clusters[[i]])
		nodes = length(unique(c(x$edgelist[clusters[[i]],1],x$edgelist[clusters[[i]],2])))
		ld[i] <- (edges-nodes+1)/(((nodes*(nodes-1))/2)-nodes+1)
		}

	names(ld) <- clusterids
	return(ld)

	}


getCommunityConnectedness <- function(x, clusterids = 1:x$numbers[3], conn = "conn", normalise = TRUE, verbose = FALSE)
	# Returns ratio of inter-community edges to intra-community edges for connectedness, or the inverse for modularity.
	# x is a "linkcomm" object.
	{
	clusters <- x$clusters[clusterids]
	intra <- NULL
	inter <- NULL
	ratio <- NULL
	if(normalise){
		avdeg <- mean(degree(x$igraph, v=V(x$igraph)))
		}
	for(i in 1:length(clusters)){
		if(verbose){
			out <- paste(c("   Completed... ",floor((i/length(clusters))*100),"%"),collapse="")
			cat(out,"\r")
			flush.console()
			}
		intra <- length(clusters[[i]])
		nodes <- unique(c(x$edgelist[clusters[[i]],1],x$edgelist[clusters[[i]],2]))
		edgeIDs <- function(y,lc){return(unique(c(which(lc$edgelist[,1]==y),which(lc$edgelist[,2]==y))))}
		alledges <- unlist(sapply(nodes,edgeIDs,lc=x))
		inter <- length(setdiff(alledges,clusters[[i]]))
		if(inter == 0){
			inter <- 1*10^-5 # Prevent infinities.
			}
		if(normalise){
			intra <- (intra)/((length(nodes)*(length(nodes)-1))/2)
			if(intra == 0){
				inter <- NULL # Trivial cluster.
			}else{
				inter <- inter/(length(nodes)*avdeg)
				}
			}
		if(conn == "conn"){
			if(length(inter)>0){
				ratio[i] <- inter/intra
			}else{
				ratio[i] <- 0
			}
		}else{
			if(length(inter)>0){
				ratio[i] <- intra/inter
			}else{
				ratio[i] <- 0
				}
			}
		}
	names(ratio) <- clusterids
	cat("\n")
	
	return(ratio)

	}


getClusterRelatedness <- function(x, clusterids = 1:x$numbers[3], hcmethod = "ward.D", cluster = TRUE, plot = TRUE, cutat = NULL, col = TRUE, pal = brewer.pal(11,"Spectral"), labels = FALSE, plotcut = TRUE, right = TRUE, verbose = TRUE, ...)
	# Returns hclust object and plots the dendrogram of cluster relatedness based on nodes.
	# Uses Jaccard coefficient to assign relatedness based on the number of shared nodes.
	# x is a "linkcomm" object.
	{
	# Check for downwards-incompatibility of ward.D:
	if(hcmethod == "ward.D" || hcmethod == "ward.D2"){
		qq <- sessionInfo()
		vers <- as.numeric(substr(paste(qq$R.version$major, qq$R.version$minor, sep="."),1,3))
		if(vers < 3.1){
			hcmethod <- "ward"
			}
		}

	nodes <- x$nodeclusters[x$nodeclusters[,2]%in%clusterids,1]
	clusters <- x$nodeclusters[x$nodeclusters[,2]%in%clusterids,2]
	numN <- length(nodes)

	emptyvec <- rep(1,(length(clusterids)*(length(clusterids)-1))/2)

	dissvec <- .C("getJaccards", as.integer(nodes), as.integer(clusters), as.integer(clusterids), as.integer(numN), dissvec = as.double(emptyvec), as.logical(verbose))$dissvec

	if(cluster){
		distmatrix <- matrix(1,length(clusterids),length(clusterids))
		distmatrix[lower.tri(distmatrix)] <- dissvec
		colnames(distmatrix) <- clusterids
		rownames(distmatrix) <- clusterids
		distobj <- as.dist(distmatrix) # Convert into 'dist' object for hclust.
		rm(distmatrix)
		cat("\n   Hierarchical clustering...\n")
		hcl <- hclust(distobj, method = hcmethod)
		hcl$order <- rev(hcl$order)
		if(length(cutat)==0){
			if(plot){
				cat("   Plotting... \n")
				plot(hcl)
				}
			return(hcl)
		}else{
			clus <- cutDendrogramAt(hcl, lc = x, cutat = cutat, plot = plot, col = col, pal = pal, labels = labels, plotcut = plotcut, right = right)

			return(clus)
			
			}
	}else{
		if(verbose){
			cat("\n")
			}
		return(dissvec)
		}

	}


cutDendrogramAt <- function(x, lc = NULL, cutat = NULL, plot = TRUE, col = TRUE, pal = brewer.pal(9,"Set1"), labels = FALSE, plotcut = TRUE, right = TRUE, verbose = TRUE, ...)
	# Returns clusters from a dendrogram after cutting at a user-chosen height.
	# x is an "hclust" object.
	{
	numM <- length(which(x$height <= cutat))

	csize <- .C("cutTreeAt", as.integer(x$merge[1:numM,1]), as.integer(x$merge[1:numM,2]), as.double(x$height[1:numM]), as.double(cutat), csize = integer(1), as.integer(numM))$csize
			
	# Read in clusters from a file.
	clus <- list()
	for(i in 1:csize){
		clus[[i]] <- scan(file = "linkcomm_metaclusters.txt", nlines = 1, skip = i-1, quiet = TRUE)
		}
	cat("\n")
	file.remove("linkcomm_metaclusters.txt")

	if(plot){
		dd <- as.dendrogram(x)
		if(col){
			crf <- colorRampPalette(pal,bias=1)
			cols <- crf(length(clus))
			cols <- sample(cols,length(clus),replace=FALSE)
			numnodes <- nrow(x$merge) + length(which(x$merge[,1]<0)) + length(which(x$merge[,2]<0))
			dd <- dendrapply(dd, .COL, height=cutat, clusters=unlist(clus), cols=cols, labels=labels, numnodes=numnodes, droptrivial = FALSE, verbose=verbose)
			assign("i",0,environment(.COL))
			assign("memb",0,environment(.COL))
			assign("first",0,environment(.COL))
			assign("left",0,environment(.COL))
			}
		if(right){dd <- rev(dd)}
		cat("\n   Plotting... \n")
		plot(dd,ylab="Height", ...)
		if(plotcut){
			abline(h=cutat,col='red',lty=2,lwd=2)
			}
		if(length(lc)>0){
			ll <- sapply(clus,length)
			maxnodes <- length(unique(lc$nodeclusters[lc$nodeclusters[,2]%in%unlist(clus[which(ll==max(ll))]),1]))
			summ <- paste("# clusters = ",length(clus),"\nLargest cluster = ",maxnodes," nodes")
			mtext(summ, line = -28)
		}else{
			summ <- paste("# clusters = ",length(clus))
			mtext(summ, side = 1)
			}
		}

	return(clus)

	}


newLinkCommsAt <- function(x, cutat = 0.5)
	# Returns a "linkcomm" object with clusters derived from cutting the dendrogram at the specified height.
	# x is a "linkcomm" object.
	{
	numM <- length(which(x$hclust$height <= cutat))

	csize <- .C("cutTreeAt", as.integer(x$hclust$merge[1:numM,1]), as.integer(x$hclust$merge[1:numM,2]), as.double(x$hclust$height[1:numM]), as.double(cutat), csize = integer(1), as.integer(numM))$csize

	if(csize == 0){
		stop("\nThere are no clusters appearing at this height; maybe try a different height.\n")
		}

	# Read in clusters from a file.
	clus <- list()
	for(i in 1:csize){
		clus[[i]] <- scan(file = "linkcomm_metaclusters.txt", nlines = 1, skip = i-1, quiet = TRUE)
		}
	cat("\n")
	file.remove("linkcomm_metaclusters.txt")

	edges <- integer.edgelist(x$edgelist)$edges

	# Extract nodes for each edge cluster.
	ecn <- data.frame()
	ee <- data.frame()
	for(i in 1:(length(clus))){
		mes<-paste(c("   Finishing up...1/3... ",floor((i/length(clus))*100),"%"),collapse="")
		cat(mes,"\r")
		flush.console()
		ee <- rbind(ee,cbind(x$edgelist[clus[[i]],],i))
		nodes <- V(x$igraph)$name[(unique(c(edges[clus[[i]],])))]
		both <- cbind(nodes,rep(i,length(nodes)))
		ecn <- rbind(ecn,both)
		}
	colnames(ecn) <- c("node","cluster")
	colnames(ee) <- c("node1","node2","cluster")

	# Extract the node-size of each edge cluster and order largest to smallest.
	ss <- NULL
	unn <- unique(ecn[,2])
	for(i in 1:length(unn)){
		mes<-paste(c("   Finishing up...2/3... ",floor((i/length(unn))*100),"%"),collapse="")
		cat(mes,"\r")
		flush.console()
		ss[i] <- length(which(ecn[,2]==unn[i]))
		}
	names(ss) <- unn
	ss <- sort(ss,decreasing=T)

	# Extract the number of edge clusters that each node belongs to.
	oo<-NULL
	unn <- unique(ecn[,1])
	for(i in 1:length(unn)){
		mes<-paste(c("   Finishing up...3/3... ",floor((i/length(unique(ecn[,1])))*100),"%"),collapse="")
		cat(mes,"\r")
		flush.console()
		oo[i]<-length(which(ecn[,1]==unn[i]))
		}
	cat("\n")
	names(oo)<-unn

	# Add nodeclusters of size 0.
	missnames <- setdiff(V(x$igraph)$name,names(oo))
	m <- rep(0,length(missnames))
	names(m) <- missnames
	oo <- append(oo,m)

	x$numbers[3] <- length(clus)
	x$nodeclusters <- ecn
	x$clusters <- clus
	x$edges <- ee
	x$pdmax <- cutat
	x$numclusters <- sort(oo,decreasing=TRUE)
	x$clustsizes <- ss
	
	return(x)

	}


getNodesIn <- function(x, clusterids = 1, type = "names")
	# x is a "linkcomm" object.
	{
	nodes <- unique(x$nodeclusters[x$nodeclusters[,2]%in%clusterids,1])
	if(type == "names"){
		return(as.character(nodes))
	}else if (type == "indices"){
		inds <- match(nodes, V(x$igraph)$name)
		return(inds)
		}
	}


getEdgesIn <- function(x, clusterids = 1, nodes = NULL, all = FALSE)
	# x is a "linkcomm" or "OCG" object.
	{
	cl <- class(x)
	switch(cl,
		linkcomm = {
	if(is.null(nodes)==TRUE){
		edges <- x$clusters[clusterids]
		return(unlist(edges))
	}else{
		if(all){ # Return edges in all of the communities to which the node(s) belong.
			clusterids <- unique(x$nodeclusters[x$nodeclusters[,1]%in%nodes,2])
			edges <- x$clusters[clusterids]
			return(unlist(edges))			
		}else{ # Return only edges directly incident upon the node(s).
			edges <- NULL
			for(i in 1:length(nodes)){
				edges <- append(edges, unique(c(which(x$edgelist[,1]==nodes), which(x$edgelist[,2]==nodes))))
				}
			return(unique(edges))
			}
		}
		},
		OCG = {
	nodes <- unique(x$nodeclusters[x$nodeclusters[,2]%in%clusterids,1])
	if(all){
		eids <- 1:nrow(x$edgelist)
		edges <- c(eids[x$edgelist[,1]%in%nodes],eids[x$edgelist[,2]%in%nodes])
		return(unique(edges))
	}else{	# Only interactions between nodes within the community.
		eids <- 1:nrow(x$edgelist)
		edges <- intersect(eids[x$edgelist[,1]%in%nodes],eids[x$edgelist[,2]%in%nodes])
		return(edges)
		}
		}
		)
	}


graph.feature <- function(x, type = "nodes", clusterids = 1:length(x$clusters), nodes = NULL, indices, features, default = 15, showall = FALSE)
	# x is a "linkcomm" object.
	# features must be same length as indices.
	# Returns named vector of node sizes or edge widths.
	{
	feat <- NULL
	if(length(nodes) > 0){
		clusterids <- unique(x$nodeclusters[x$nodeclusters[,1]%in%nodes,2])
		}
	clusters <- x$clusters[clusterids]
	if(showall){
		# Add single edge "clusters".
		single <- setdiff(1:x$numbers[1],unlist(clusters))
		ll <- length(clusters)
		for(i in 1:length(single)){
			clusters[[(i+ll)]] <- single[i]
			}
		}
	if(type == "edges"){
		edges <- x$edgelist[unlist(clusters),]
		ig <- graph.edgelist(edges, directed=x$directed)
		if(length(unlist(clusters)) >= nrow(x$edgelist)){ # No changes in index order needs to be made.
			feat <- rep(default, length(E(ig)))
			inds <- match((indices-1), E(ig))
			feat[inds] <- features
			names(feat) <- E(ig)
		}else{
			feat <- rep(default, length(E(ig)))
			inds <- match(indices, unlist(clusters)) # Where are the edge indices in the re-ordered sub-network.
			feat[inds] <- features
			names(feat) <- E(ig)
			}
	}else if(type == "nodes"){
		nodes <- V(x$igraph)$name[indices]
		edges <- x$edgelist[unlist(clusters),]
		ig <- graph.edgelist(edges, directed=x$directed)
		if(length(unlist(clusters)) >= nrow(x$edgelist)){ # No changes in index order needs to be made.
			feat <- rep(default, length(V(ig)))
			feat[indices] <- features
			names(feat) <- V(ig)$name
			feat <- feat[match(names(feat), V(ig)$name)]
		}else{
			feat <- rep(default, length(V(ig)))
			inds <- match(nodes, V(ig)$name)
			feat[inds] <- features
			names(feat) <- V(ig)$name
			feat <- feat[match(names(feat), V(ig)$name)]
			}
		}

	return(feat)

	}


get.community.overlaps <- function(x){
	# x is a linkcomm or OCG object.
	tt <- list()
	for(i in 1:x$numbers[3]){tt[[i]]<-NA}
	for(i in 1:x$numbers[2]){
		ww <- which.communities(x, nodes = V(x$igraph)$name[i])
		if(length(ww)>0){
			for(j in 1:length(ww)){
				tt[[ww[j]]] <- append(tt[[ww[j]]], ww)
				tt[[ww[j]]] <- tt[[ww[j]]][-which(tt[[ww[j]]] == ww[j])]
				}
			}
		}
	# Clean up NA's and make communities unique.
	for(i in 1:length(tt)){
		if(length(which(tt[[i]]>0))>0){
			tt[[i]] <- tt[[i]][-which(is.na(tt[[i]]))]
			tt[[i]] <- sort(unique(tt[[i]]))
			}
		}
	return(tt)
	}


get.shared.nodes <- function(x, comms){
	# x is a linkcomm or OCG object.
	# comms is a vector of community IDs.
	qq <- list()
	for(i in 1:length(comms)){
		qq[[i]] <- getNodesIn(x, clusterids=comms[i])
		}
	nn <- Reduce(intersect, qq)
	return(nn)
	}


meta.communities <- function(x, hcmethod = "ward.D", deepSplit = FALSE)
	{
	# x is a linkcomm or OCG object.

	dissvec <- getClusterRelatedness(x, cluster = FALSE)

	distmatrix <- matrix(1,x$numbers[3], x$numbers[3])
	distmatrix[lower.tri(distmatrix)] <- dissvec
	colnames(distmatrix) <- 1:x$numbers[3]
	rownames(distmatrix) <- 1:x$numbers[3]
	cat("   Hierarchical clustering...\n")

	# Check for downwards-incompatibility of ward.D:
	if(hcmethod == "ward.D" || hcmethod == "ward.D2"){
		qq <- sessionInfo()
		vers <- as.numeric(substr(paste(qq$R.version$major, qq$R.version$minor, sep="."),1,3))
		if(vers < 3.1){
			hcmethod <- "ward"
			}
		}

	hcl <- hclust(as.dist(distmatrix), method = hcmethod)

	scl <- cutreeHybrid(hcl, distM = distmatrix, deepSplit = deepSplit)

	scl <- scl$labels
	names(scl) <- 1:x$numbers[3]
	# Zero indicates community did not get clustered.

	# Modify original linkcomm object to reflect new communities.
	nc <- length(unique(scl))
	x$nodeclusters[,2] <- scl[match(x$nodeclusters[,2], names(scl))]
	# Remove duplicates.
	dd <- duplicated(x$nodeclusters)
	wd <- which(dd==TRUE)
	if(length(wd) > 0){
		x$nodeclusters <- x$nodeclusters[-wd,]
		}
	
	# Number of nodes in each cluster.
	nn <- NULL
	for(i in 1:nc){
		nn <- append(nn, length(which(x$nodeclusters[,2]==i)))
		}
	x$clustsizes <- nn
	names(x$clustsizes) <- 1:nc
	# Need a new edge cluster list which merges edges into new clusters.
	if(class(x)=="linkcomm"){
		ecl <- list()
		for(i in 1:nc){
			# Get old communities that belong to this new cluster.
			oc <- as.integer(names(scl)[which(scl==i)])
			ecl[[i]] <- sort(unlist(x$clusters[oc]))
			}
		x$clusters <- ecl
		x$edges[,3] <- scl[match(x$edges[,3], names(scl))]
		}

	# Extract the number of edge clusters that each node belongs to.
	unn <- unique(x$nodeclusters[,1])

	iecn <- as.integer(as.factor(x$nodeclusters[,1]))
	iunn <- unique(iecn)
	lunn <- length(iunn)
	nrows <- nrow(x$nodeclusters)

	oo <- rep(0,lunn)

	verbose <- TRUE

	oo <- .C("getNumClusters", as.integer(iunn), as.integer(iecn), counts = as.integer(oo), as.integer(lunn), as.integer(nrows), as.logical(verbose))$counts
	cat("\n")

	names(oo) <- unn
	
	x$numclusters <- oo
	# Finally, change the number of communities.
	x$numbers[3] <- nc

	return(x)

	}









