igraphKamadaKawai <- function(n)
{
	if(length(n) != 1 || mode(n) != "numeric") stop("Input n must be a single number")
	g <- graph.tree(n)
	g <- g + vertex(1)
	g <- g + edge(n+1, 1)
	locations_this_network <- layout.kamada.kawai(graph = g)

	return(list(locations = locations_this_network, graph = g, initialPoint = n+1))
}
iterativeTreeLayout <- function(n)
{
	# minimum distance entre 2 segments (version for quick programming (quick and dirty)
	ds <- function(x1,x2,y1,y2,x3,x4,y3,y4)
	{
		xp <- seq(x1,x2,length=10)
		yp <- seq(y1,y2,length=10)
		xq <- seq(x3,x4,length=10)
		yq <- seq(y3,y4,length=10)
		dps <- distpq(xp,yp,xq,yq)
		mds <- min(dps)
		return(mds)
	}

	distpq <- function(xp,yp,xq,yq)
	{
		np <- length(xp)
		nq <- length(xq)
		dist <- matrix(0, np, nq)  
		dist <- dist + outer(xp,xq,"-")^2
		dist <- dist + outer(yp,yq,"-")^2
		dist <- sqrt(dist)
		return(dist)
	}

	tree <- matrix(NA,n,7)
	colnames(tree) <- c("indseg","x","y","ancestor","nbtrib","length","theta")
	tree[1,1] <- 1 
	tree[1,2] <- 0
	tree[1,3] <- 1
	tree[1,4] <- 0
	tree[1,5] <- 0
	tree[1,6] <- 1
	tree[1,7] <- 0

	i <- 2
	tree[i,1] <- 2  # i
	tree[i,4] <- 1
	tree[i,5] <- 0
	tree[tree[i,4],5] <- 1
	length <- runif(1,0.9,1.1)
	theta <- runif(1,-pi/4,+pi/4)                   
	tree[i,2] <- tree[tree[i,4],2] + length * sin(theta)
	tree[i,3] <- tree[tree[i,4],3] + length * cos(theta)           
	tree[i,6] <- length
	tree[i,7] <- theta              
	  
	for (i in 3:n) 
	{ 
		tree[i,4] <- sample(seq(1,(i-1),1),1)
		while(tree[tree[i,4],5]>=2) {tree[i,4] <- sample(seq(1,(i-1),1),1)} 

		tree[i,5] <- 0
		# tree[tree[i,4],5] <- tree[tree[i,4],5] + 1  
		theta <- runif(1,-pi/4,+pi/4)
		length <- runif(1,0.9,1.1)
		tree[i,2] <- tree[tree[i,4],2] + length * sin(theta)
		tree[i,3] <- tree[tree[i,4],3] + length * cos(theta)           
		tree[i,6] <- length
		tree[i,7] <- theta

		# calcul de la distance minimum avec les autres segments except
		# son "pere" et son eventuel "frere" si il existe  
		dsi <- Inf
		for (k in 2:(i-1))
		{
			if ( k != tree[i,4] & tree[k,4] != tree[i,4] ) 
			{
				dsi <- min(dsi,ds(tree[i,2],tree[tree[i,4],2],
					tree[i,3],tree[tree[i,4],3],
					tree[k,2],tree[tree[k,4],2],
					tree[k,3],tree[tree[k,4],3]))
			} 
		}
		crit <- dsi 

		nbloops <- 1   
		while(crit <= 0.2)
		{
			tree[i,4] <- sample(seq(1,(i-1),1),1)
			while(tree[tree[i,4],5]>=2)
			{
				tree[i,4] <- sample(seq(1,(i-1),1),1)
			} 

			tree[i,5] <- 0
			# tree[tree[i,4],5] <- tree[tree[i,4],5] + 1  
			theta <- runif(1,-pi/4,+pi/4)
			length <- runif(1,0.9,1.1)
			tree[i,2] <- tree[tree[i,4],2] + length * sin(theta)
			tree[i,3] <- tree[tree[i,4],3] + length * cos(theta)           
			tree[i,6] <- length
			tree[i,7] <- theta

			dsi <- Inf
			for (k in 2:(i-1))
			{
				if ( k != tree[i,4] & tree[k,4] != tree[i,4] ) 
				{
					dsi <- min(dsi,ds(tree[i,2],tree[tree[i,4],2],
							  tree[i,3],tree[tree[i,4],3],
							  tree[k,2],tree[tree[k,4],2],
							  tree[k,3],tree[tree[k,4],3]))
				} 
			}
			crit <- dsi
			nbloops <- nbloops + 1
			if(nbloops > 1000) 
			{ 
				stop("no more possibilities to add a new branch")
			}
		}  
		tree[tree[i,4],5] <- tree[tree[i,4],5] + 1
	}
	edgelist <- cbind(tree[,"ancestor"]+1, 2:(nrow(tree)+1))
	return(list(graph=graph.edgelist(edgelist, directed=FALSE), locations = rbind(c(0,0), cbind(tree[,"x"], tree[, "y"])), initialPoint = 1))
}

createSSN <- function(n, obsDesign, predDesign = noPoints, path, importToR = FALSE, treeFunction = igraphKamadaKawai)
{
#	if(!require(igraph))
#	{
#		stop("simulate requires access to the igraph package")
#	}
	if(missing(obsDesign)) stop("Input obsDesign cannot be missing")
	if(missing(path)) stop("Path cannot be missing")
	if(missing(n))
	{
		stop("Input n cannot be missing")
	}
	if(length(path) != 1) stop("Please enter a single path")
	info <- file.info(path)
	isdir <- info$isdir
	if(is.na(isdir))
	{
		dir.create(path)
	}
	else if(isdir == FALSE)
	{
		stop("Unable to create directory")
	}
	
	old_wd <- getwd()
	on.exit(setwd(old_wd))
	setwd(path)

	#construct graph, get edges
	
	n_networks <- length(n)
	edges <- vector(mode="list", length=n_networks)
	tree.graphs <- edges
	locations <- edges
	rids <- edges
	initial_points <- n_edges <- vector(mode="numeric", length = n_networks)
	edge_lengths <- edges
	cumulative_nedges <- 0
	#get the maximum and minimum x values, so we can offset the networks and have them spread out in the x direction
	max_x <- vector(mode="numeric", length=n_networks)
	min_x <- vector(mode="numeric", length=n_networks)
	for(i in 1:n_networks) 
	{
		graph <- treeFunction(n[i])

		edges_this_network <- get.edgelist(graph$graph, names=FALSE)
		reordering <- order(edges_this_network[,1])
		
		edges_this_network <- edges_this_network[reordering,]
		locations_this_network <- graph$locations[reordering]
		tree.graphs[[i]] <- graph.edgelist(edges_this_network)
		
		initial_points[i] <- graph$initialPoint
		locations_this_network <- graph$locations
		locations[[i]] <- locations_this_network
		
		edges[[i]] <- edges_this_network
		
		
		n_edges[i] = nrow(edges_this_network)
		
		#calculate edge lengths
		edge_lengths_function <- function(indicies)
		{
			return(sqrt(sum((locations_this_network[indicies[1],] - locations_this_network[indicies[2],])^2)))
		}
		edge_lengths[[i]] <- apply(edges_this_network, 1, edge_lengths_function)
		rids[[i]] <- (1:n_edges[i]) + cumulative_nedges
		names(edge_lengths[[i]]) <- rids[[i]]
		cumulative_nedges <- cumulative_nedges + n_edges[i]
		min_x[i] <- min(locations_this_network)
		max_x[i] <- max(locations_this_network)
	}
	#do the offsetting
	cumulative_x <- max_x[1]
	if(n_networks > 1)
	{
		for(i in 2:n_networks)
		{
			locations[[i]][,1] <- locations[[i]][,1] + cumulative_x - min_x[i] + 0.1
			cumulative_x <- cumulative_x + max_x[i] - min_x[i] + 0.1
		}
	}
	#set up storage for the edges data to be written by writeSpatialShape
	spatial_edges <- vector(mode="list", length=sum(n_edges))
	cumulative_nedges <- 0
	for(netid in 1:n_networks)
	{
		locations_this_network <- locations[[netid]]
		edges_this_network <- edges[[netid]]
		for(edge.index in 1:n_edges[netid])
		{
			edge <- edges_this_network[edge.index,]
			first.location = locations_this_network[edge[1],]
			second.location = locations_this_network[edge[2],]
			#Reverse the locations, because we're dealing with them branching OUT from the root point, whereas we actually want the opposite directions when we save them
			spatial_edges[[edge.index + cumulative_nedges]] = Lines(list(Line(rbind(second.location, first.location))), ID=as.character(edge.index+cumulative_nedges))
		}
		cumulative_nedges = cumulative_nedges + n_edges[netid]
	}
	sl <- SpatialLines(spatial_edges)
	
	#compute the upDist stuff
	edge_updist <- vector(mode="list", length=n_networks)
	line_data <- data.frame()
	shreve <- vector(mode="list", length=n_networks)
	for(netid in 1:n_networks)
	{
		rids_this_network <- rids[[netid]]
		edges_this_network <- edges[[netid]]
		#calculate upstream distances
		edge_updist_this_network = 0
		#point numbers
		known_points = initial_points[netid]
		remaining_edges <- edges_this_network
		edge_lengths_this_network <- edge_lengths[[netid]]
		remaining_edge_lengths <- edge_lengths_this_network
		known_rids <- c()
		remaining_rids <- rids_this_network
		while(TRUE)
		{
			can_calculate <- (remaining_edges[,1] %in% known_points) & (!(remaining_edges[,2] %in% known_points))
			upstream_point_indicies <- remaining_edges[can_calculate, 2]
			downstream_point_indicies <- remaining_edges[can_calculate, 1]
			
			edge_updist_this_network <- c(edge_updist_this_network, edge_updist_this_network[match(downstream_point_indicies, known_points)] + remaining_edge_lengths[can_calculate])
			known_points <- c(known_points, upstream_point_indicies)
			
			remaining_edges <- remaining_edges[!can_calculate,,drop=FALSE]
			remaining_edge_lengths <- remaining_edge_lengths[!can_calculate]

			known_rids <- c(known_rids, remaining_rids[can_calculate])
			remaining_rids <- remaining_rids[!can_calculate]

			if(length(remaining_edges) == 0) break
		}
		#Remove the dummy value 0 which initiated the computation. And reorder because we're going to throw away known_rids
		edge_updist_this_network <- edge_updist_this_network[-1][order(known_rids)]
		names(edge_updist_this_network) <- sort(known_rids)
		edge_updist[[netid]] <- edge_updist_this_network
		
		#compute shreves stream order for each segment
		shreve_this_network <- vector(mode="numeric", length=n_edges[netid])
		#We start off with the terminal points having known shreve values
		is_initial <- !(edges_this_network[,2] %in% edges_this_network[,1])
		
		known_points <- edges_this_network[which(is_initial),2]
		remaining_points <- edges_this_network[which(!is_initial),2]
		shreve_values_points <- rep(1, length(known_points))
		
		shreve_values_rids <- c()
		known_rids <- c()
		#Here we want to pretend that an rid is the index of an edge. But multiple networks add an offset. So subtract it for this loop. 
		remaining_rids <- rids_this_network - min(rids_this_network) + 1
		
		remaining_edges <- edges_this_network
		while(TRUE)
		{
			#First edges
			can_calculate <- (remaining_edges[,2] %in% known_points)

			shreve_values_rids <- c(shreve_values_rids, shreve_values_points[match(remaining_edges[can_calculate, 2], known_points)])
			remaining_edges <- remaining_edges[!can_calculate, , drop =FALSE]
			known_rids <- c(known_rids, remaining_rids[can_calculate])
			remaining_rids <- remaining_rids[!can_calculate]
			
			if(length(remaining_edges) == 0) break
			#then points
			can_calculate <- !(remaining_points %in% remaining_edges[,1])
			
			#are the 1/2 upstream edges already calculated?
			can_calculate_function <- function(index)
			{
				#Which edges start at this point?
				relevant_edges <- which(edges_this_network[,1] == remaining_points[index])
				#Are they already in the calculated set?
				return(all(relevant_edges %in% known_rids))
			}
			#this means that some of the true values can turn false
			can_calculate[can_calculate] <- sapply(which(can_calculate), can_calculate_function)
			
			calculate_shreve <- function(index)
			{
				#Which edges start at this point?
				relevant_edges <- which(edges_this_network[,1] == remaining_points[index])
				#Are they already in the calculated set?
				relevant_known_rids <- match(relevant_edges, known_rids)
				return(sum(shreve_values_rids[relevant_known_rids]))
			}
			new_shreve_values_points <- sapply(which(can_calculate), calculate_shreve)
			
			shreve_values_points <- c(shreve_values_points, new_shreve_values_points)
			known_points <- c(known_points, remaining_points[can_calculate])
			remaining_points <- remaining_points[!can_calculate]
		}
		#Reorder, so we're back to a standard ordering and match up with the edge_updist_this_network calculations
		known_rids_ordering <- order(known_rids)
		shreve_values_rids <- shreve_values_rids[known_rids_ordering]
		
		#Now work out the additive function values. 
		additive_function_values <- shreve_values_rids / max(shreve_values_rids)
		#for(point_index in 0:n[netid])
		#{
			#output_edges = which(point_index == edges_this_network[,1])
			#if(length(output_edges) == 0) next
			#additive_function_values[output_edges] <- additive_function_values[output_edges] / sum(additive_function_values[output_edges])
		#}
		#Abbreviate name for additive_function_value because of 10 character limit for column names (or so it seems)
		additional_line_data <- data.frame(rid = rids_this_network, netID = as.factor(rep(netid, n_edges[netid])), upDist = edge_updist_this_network, shreve = shreve_values_rids, Length=edge_lengths_this_network, addfunccol = additive_function_values)
		rownames(additional_line_data) <- rids_this_network
		line_data <- rbind(additional_line_data, line_data)
	}
	sldf <- SpatialLinesDataFrame(sl, data=line_data, match.ID = TRUE)
	writeSpatialShape(sldf, "edges")

	#Change the working directory back in case the path was specified as relative. We only needed it for the writeSpatialShapes calls
	setwd(old_wd)
	
	#create the binary ID files, in the form <net-name>.dat, e.g.
	##################
	#"rid", "binaryID"
	#362,1
	#368,10
	#378,11
	##################
	binary_ids_tables <- list()
	for(netid in 1:n_networks)
	{
		#The edges for which we don't yet have a binary ID
		remaining_edges <- edges[[netid]]
		
		remaining_points <- unique(remaining_edges)
		#known_point_indicies gives the indicies of the points for which we know the binary ID. We're twisting things a little here because we're really interested in edges, not points. But for a tree every edge can be treated as its end point, and it's more convenient. 
		known_point_indicies <- c(initial_points[netid])
		#apparently binary IDS
		known_binaryids <- c("")
		known_rids <- c()
		remaining_rids <- rids[[netid]]
		while(TRUE)
		{
			next_points <- (remaining_edges[,2] %in% remaining_points) & (remaining_edges[,1] %in% known_point_indicies)
			upstream_points <- remaining_edges[next_points, 2]
			downstream_points <- remaining_edges[next_points, 1]
			
			additional_binary_bits <- c()
			counter <- 1
			while(counter <= length(downstream_points))
			{
				if(counter == length(downstream_points)	|| downstream_points[counter] != downstream_points[counter+1])
				{
					additional_binary_bits <- c(additional_binary_bits, "1")
					counter <- counter + 1
				}
				else 
				{
					additional_binary_bits <- c(additional_binary_bits, "0", "1")
					counter <- counter + 2
				}
			}
			
			binaryid_indicies <- match(downstream_points, known_point_indicies)
			previous_binary_ids <- known_binaryids[binaryid_indicies]
			
			#This works because edges was initially sorted by downstream point.
			known_binaryids <- c(known_binaryids, paste(previous_binary_ids, additional_binary_bits, sep=""))
			known_point_indicies <- c(known_point_indicies, upstream_points)
			known_rids <- c(known_rids, remaining_rids[next_points])
			remaining_edges <- remaining_edges[!next_points,, drop=FALSE]
			remaining_rids <- remaining_rids[!next_points]
			remaining_points <- remaining_points[remaining_points %in% remaining_edges[,2]]
			
			if(length(remaining_rids) == 0) break
		}
		#strip off the extraneous initial value
		known_binaryids <- known_binaryids[-1]
		binary_ids <- known_binaryids[order(known_rids)]
		binary_ids_tables[[netid]] <- data.frame(rid = rids[[netid]], binaryID = binary_ids)
		if(length(unique(binary_ids_tables[[netid]]$binaryID)) != length(binary_ids_tables[[netid]]$binaryID)) stop("Internal error")
		write.table(binary_ids_tables[[netid]], file = file.path(path, paste("netID", netid, ".dat", sep="")), col.names=T, sep=",", row.names=FALSE)
	}
	setwd(path)
	
	#Come up with the distance matrices between the vertices of the networks. This might be needed to generate the observation / design points. This starts out being organised by EDGES, rather than vertices. That is, matrix[i, j] gives the distance between the two down-stream vertices of edges i and j. But once it's been created and populated we change over to organising it by vertices - 
	distance_matrices <- list()
	for(netid in 1:n_networks)
	{
		edge_updist_this_network <- edge_updist[[netid]]
		edges_this_network <- edges[[netid]]
		binary_id_table <- binary_ids_tables[[netid]]
		
		distance_matrix <- matrix(0, nrow(binary_id_table), nrow(binary_id_table))
		colnames(distance_matrix) <- rownames(distance_matrix) <- binary_id_table$rid
		partial_match_function <- function(binary_id1, binary_id2)
		{
			min_len <- min(nchar(binary_id1), nchar(binary_id2))
			for(j in 1:min_len)
			{
				if(substr(binary_id1, j, j) != substr(binary_id2, j, j)) return(j-1)
			}
			return(min_len)
		}
		character_binary_ids <- as.character(binary_id_table$binaryID)
		for(i in 1:nrow(binary_id_table))
		{
			current_binary_id <- binary_id_table$binaryID[i]
			current_rid <- as.character(binary_id_table$rid[i])
			current_updist <- edge_updist_this_network[current_rid]
			
			matching_characters <- sapply(character_binary_ids, partial_match_function, as.character(current_binary_id))
			matching_substring <- substr(binary_id_table$binaryID, 1, matching_characters)
			#The rid values (key value for an edge) giving the earliest common point in the tree structure between the specific current_binary_id and all other binary_id values
			indices <- match(matching_substring, binary_id_table$binaryID)
			if(any(is.na(indices))) stop("Internal Error")
			downstream_rids <- binary_id_table$rid[indices]
			downstream_updists <- edge_updist_this_network[as.character(downstream_rids)]
			
			distance_matrix[as.character(current_rid),] <- pmax(current_updist - downstream_updists, rep(0, length(downstream_updists)))
		}
		reindex <- match(binary_id_table$rid, rids[[netid]])
		distance_matrix <- distance_matrix[reindex, reindex]
		distance_matrices[[netid]] <- distance_matrix + t(distance_matrix)
	}
	
	#come up with the prediction sites
	obs_sites <- obsDesign(tree.graphs, edge_lengths, locations, edge_updist, distance_matrices)
	#and now the observed sites
	pred_sites <- predDesign(tree.graphs, edge_lengths, locations, edge_updist, distance_matrices)

	max_observed_locID <- max(unlist(lapply(obs_sites, function(x) max(x$locID))))
	for(i in 1:length(pred_sites)) 
	{
		pred_sites[[i]]$locID <- pred_sites[[i]]$locID + max_observed_locID
	}
	
	n_obs_sites <- unlist(lapply(obs_sites, function(x) return(dim(x)[1])))
	n_pred_sites <- unlist(lapply(pred_sites, function(x) return(dim(x)[1])))
	
	
	#We need the upstream distance stuff for the sites too
	sites_data <- data.frame()
	combined_site_location_data <- c()
	
	pred_data <- data.frame()
	combined_pred_location_data <- c()
	
	cumulative_pids <- 0
	for(netid in 1:n_networks)
	{
		edge_lengths_this_network <- edge_lengths[[netid]]
		edges_this_network <- edges[[netid]]
		edge_updist_this_network <- edge_updist[[netid]]
		locations_this_network <- locations[[netid]]
		n_locations_this_network <- n_obs_sites[netid] + n_pred_sites[netid]
		rids_this_network <- rids[[netid]]

		pred_sites_this_network <- pred_sites[[netid]]
		obs_sites_this_network <- obs_sites[[netid]]
		
		f <- function(row)
		{
			rid <- as.character(row[1])
			edge_id <- match(rid, rids_this_network)
			proportion <- as.numeric(row[2])
			downstream_point <- edges_this_network[edge_id, 1]
			upstream_point <- edges_this_network[edge_id, 2]
			downstream_location <- locations_this_network[downstream_point, ]
			upstream_location <- locations_this_network[upstream_point, ]
			location <- downstream_location + proportion * (upstream_location - downstream_location)
			ret <- c(location, (sqrt(sum((location - downstream_location)^2)) + edge_updist_this_network[rid] - edge_lengths_this_network[rid]))
			names(ret) <- c("NEAR_X", "NEAR_Y", "upDist")
			return(ret)
		}
		obs_location_data <- data.frame(rid=obs_sites_this_network$edge, ratio=obs_sites_this_network$ratio, locID = obs_sites_this_network$locID, stringsAsFactors=FALSE)
		pred_location_data <- data.frame(rid=pred_sites_this_network$edge, ratio=pred_sites_this_network$ratio, locID = pred_sites_this_network$locID, stringsAsFactors=FALSE)
		
		#aggregate some more informative data into location_data_this_network
		if(n_locations_this_network > 0)
		{

			pred_location_data_this_network <- t(apply(pred_location_data, 1, f))
			if(length(pred_location_data_this_network) > 0){
                        	colnames(pred_location_data_this_network) <- c("NEAR_X", "NEAR_Y", "upDist")
			} else {
				pred_location_data_this_network <- matrix(0, 0, 3)
				colnames(pred_location_data_this_network) <- c("NEAR_X", "NEAR_Y", "upDist")
			}
			obs_location_data_this_network <- t(apply(obs_location_data, 1, f))
			if(length(obs_location_data_this_network) > 0){
                        	colnames(obs_location_data_this_network) <- c("NEAR_X", "NEAR_Y", "upDist")
			} else	{
				obs_location_data_this_network <- matrix(0, 0, 3)
				colnames(obs_location_data_this_network) <- c("NEAR_X", "NEAR_Y", "upDist")
			}
			if(n_obs_sites[netid] > 0)
			{
				obs_pids <- (1:n_obs_sites[netid]) + cumulative_pids
			}
			else obs_pids <- integer(0)
			if(n_pred_sites[netid] > 0)
			{
				pred_pids <- n_obs_sites[netid] + (1:n_pred_sites[netid]) + cumulative_pids
			}
			else pred_pids <- integer(0)
		}
		else
		{
			obs_location_data_this_network <- pred_location_data_this_network <- data.frame(NEAR_X=numeric(0), NEAR_Y = numeric(0), upDist = numeric(0))
			obs_pids <- integer(0)
			pred_pids <- integer(0)
		}
		
		cumulative_pids <- cumulative_pids + n_locations_this_network
	
		#Specify the data matrix		
		#Careful with line_data[site_rids,"shreve"], the site_rids are numbers but used to look up row NAMES, not actual rows. So line_data[1, "shreve"] is not the first row. 
		obs_data_this_network <- data.frame(locID = obs_location_data[,"locID"], upDist = obs_location_data_this_network[,"upDist"], pid=obs_pids, netID= rep(netid, length(obs_pids)), rid=obs_location_data[,"rid"], ratio=obs_location_data[,"ratio"], shreve = line_data[match(obs_location_data[,"rid"], line_data[,"rid"]),"shreve"], addfunccol = line_data[match(obs_location_data[,"rid"], line_data[,"rid"]),"addfunccol"], stringsAsFactors=FALSE)
		#Extra covariate data?
		if(ncol(obs_sites_this_network) > 3)
		{
			obs_data_this_network <- cbind(obs_data_this_network, obs_sites_this_network[,-match(c("edge", "ratio", "locID"), colnames(obs_sites_this_network)), drop=FALSE])
		}
		rownames(obs_data_this_network) <- obs_pids
		rownames(obs_location_data_this_network) <- obs_pids
		
		pred_data_this_network <- data.frame(locID = pred_location_data[,"locID"], upDist = pred_location_data_this_network[,"upDist"], pid=pred_pids, netID= rep(netid, length(pred_pids)), rid=pred_location_data[,"rid"], ratio=pred_location_data[,"ratio"], shreve = line_data[match(pred_location_data[,"rid"], line_data[,"rid"]),"shreve"], addfunccol = line_data[match(pred_location_data[,"rid"], line_data[,"rid"]),"addfunccol"], stringsAsFactors=FALSE)
		if(ncol(pred_sites_this_network) > 3)
		{
			pred_data_this_network <- cbind(pred_data_this_network, pred_sites_this_network[,-match(c("edge", "ratio", "locID"), colnames(pred_sites_this_network)), drop=FALSE])
		}
		rownames(pred_data_this_network) <- pred_pids
		rownames(pred_location_data_this_network) <- pred_pids
		
		if(n_obs_sites[netid] > 0)
		{
			sites_data <- rbind(obs_data_this_network[1:n_obs_sites[netid],,drop=FALSE], sites_data)
			combined_site_location_data <- rbind(obs_location_data_this_network[1:n_obs_sites[netid],,drop=FALSE], combined_site_location_data)
		}
		
		if(n_pred_sites[netid] > 0)
		{
			pred_data <- rbind(pred_data_this_network[(1:n_pred_sites[netid]),,drop=FALSE], pred_data)
			combined_pred_location_data <- rbind(pred_location_data_this_network[(1:n_pred_sites[netid]),,drop=FALSE],			combined_pred_location_data)
		}
	}
	#save the actual points, and attach sites_data as the related data
	if(length(combined_site_location_data) == 0) stop("At least one observation site must be present")
	sites <- SpatialPointsDataFrame(combined_site_location_data[,c("NEAR_X", "NEAR_Y"),drop=FALSE], sites_data, match.ID = TRUE)
	writeSpatialShape(sites, "sites")
	
	#similarly for the pred points
	if(length(combined_pred_location_data) > 0)
	{
		preds <- SpatialPointsDataFrame(combined_pred_location_data[,c("NEAR_X", "NEAR_Y"),drop=FALSE], pred_data, match.ID = TRUE)
		writeSpatialShape(preds, "preds")
	}
	

	#Change the working directory back in case the path was specified as relative. We only needed it for the writeSpatialShapes calls
	setwd(old_wd)
	if(importToR)
	{
		if(sum(n_pred_sites) > 0) return(importSSN(path, predpts="preds", o.write=TRUE))
		else return(importSSN(path, o.write=TRUE))
	}
	else return(invisible(NULL))
}
