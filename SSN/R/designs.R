replication.function <- function(design.function, replications = 1, rep.variable = "Time", rep.values)
{
		replicated.function <- function(tree.graphs, edge_lengths, locations, edge_updist, distance_matrices)
		{
			n_networks <- length(tree.graphs)
			unreplicated <- design.function(tree.graphs, edge_lengths, locations, edge_updist, distance_matrices)
			replicated <- unreplicated
			if(replications > 1)
			{
				for(netid in 1:n_networks)
				{
					for(i in 1:(replications-1))
					{
						replicated[[netid]] <- rbind(replicated[[netid]], unreplicated[[netid]])
					}
					replicated[[netid]][,rep.variable] <- rep(rep.values, each=nrow(unreplicated[[netid]]))
					replicated[[netid]]$locID <- rep(unreplicated[[netid]]$locID, times=replications)
				}
			}
			return(replicated)
		}
		return(replicated.function)
}
systematicDesign <- function(spacing, replications=1, rep.variable = "Time", rep.values)
{
	if(missing(rep.values)) rep.values <- 1:replications
	if(replications != length(rep.values))
	{
		stop("Input rep.values must contain one element for each replication")
	}
	design.function <- function(tree.graphs, edge_lengths, locations, edge_updist, distance_matrices)
	{
		if(length(spacing) == 1) spacing <- rep(spacing, length(tree.graphs))
		
		if(length(spacing) != length(tree.graphs))
		{
			stop("Dimension mismatch: Input spacing must contain one number, or one number for each network")
		}
		n_networks <- length(tree.graphs)
		result <- vector(mode="list", length=length(n_networks))
		cumulative_locID <- 0
		for(netid in 1:n_networks)
		{
			spacing_this_network <- spacing[netid]
			
			graph <- tree.graphs[[netid]]
			edge_lengths_this_network <- edge_lengths[[netid]]
			rids <- names(edge_updist[[netid]])
			edges_this_network <- get.edgelist(graph)
			points_this_network <- sort(unique(as.numeric(edges_this_network)))
			
			positions_per_segment <- vector(mode="list", length=length(edge_lengths_this_network))
			done_points <- !(points_this_network %in% edges_this_network[,2])
			done_segments <- c()
			segment_remaining <- c()
			while(length(done_segments) != nrow(edges_this_network))
			{
				can_calculate <- done_points[match(edges_this_network[,1], points_this_network)]
				can_calculate[done_segments] <- FALSE
				can_calculate_indices <- which(can_calculate)
				if(!any(can_calculate)) stop("Internal error")
				for(index in can_calculate_indices)
				{
					edge <- edges_this_network[index,]
					remaining <- segment_remaining[match(match(edge[1], edges_this_network[,2]), done_segments)]
					if(is.null(remaining)) remaining <- spacing_this_network
					edge_length <- edge_lengths_this_network[index]
					if(edge_length + remaining < spacing_this_network)
					{
						segment_remaining <- c(segment_remaining, edge_length + remaining)
					}
					else
					{
						positions_per_segment[[index]] <- seq(spacing_this_network - remaining, edge_length, by=spacing_this_network)
						segment_remaining <- c(segment_remaining, edge_length - max(positions_per_segment[[index]]))
					}
					done_segments <- c(done_segments, index)
					done_points[match(edge[2], points_this_network)] <- TRUE
				}
			}
			proportions_per_segment <- positions_per_segment
			for(i in 1:length(proportions_per_segment)) proportions_per_segment[[i]] <- proportions_per_segment[[i]] / edge_lengths_this_network[i]
			
			unreplicated <- data.frame(edge = rep(rids, times=unlist(lapply(proportions_per_segment, length))), ratio = unlist(proportions_per_segment), stringsAsFactors=FALSE)
			unreplicated$locID <- 1:nrow(unreplicated) + cumulative_locID
			cumulative_locID <- cumulative_locID + nrow(unreplicated)
			
			result[[netid]] <- unreplicated
		}
		return(result)
	}
	return(replication.function(design.function, replications, rep.variable, rep.values))
}
binomialDesign <- function(n, replications=1, rep.variable = "Time", rep.values)
{
	if(missing(rep.values)) rep.values <- 1:replications
	if(replications != length(rep.values))
	{
		stop("Input rep.values must contain one element for each replication")
	}
	design.function <- function(tree.graphs, edge_lengths, locations, edge_updist, distance_matrices)
	{
		if(length(tree.graphs) != length(n))
		{
			stop("Dimension mismatch: Input n must contain one number for each network")
		}
		n_networks <- length(tree.graphs)
		result <- vector(mode="list", length=length(n_networks))
		cumulative_locID <- 0
		for(netid in 1:n_networks)
		{
			indices <- sample(1:length(edge_lengths[[netid]]), prob = edge_lengths[[netid]] / sum(edge_lengths[[netid]]), replace=TRUE, size=n[netid])
			rids <- names(edge_updist[[netid]])
			proportion <- runif(n[netid])
			
			unreplicated <- data.frame(edge = rids[indices], ratio = proportion, stringsAsFactors=FALSE)
			unreplicated$locID <- 1:nrow(unreplicated) + cumulative_locID
			cumulative_locID <- cumulative_locID + nrow(unreplicated)
			
			result[[netid]] <- unreplicated
		}
		return(result)
	}
	return(replication.function(design.function, replications, rep.variable, rep.values))
}
poissonDesign <- function(lambda, replications=1, rep.variable = "Time", rep.values)
{
	if(missing(rep.values)) rep.values <- 1:replications
	if(replications != length(rep.values))
	{
		stop("Input rep.values must contain one element for each replication")
	}
	design.function <- function(tree.graphs, edge_lengths, locations, edge_updist, distance_matrices)
	{
		if(length(tree.graphs) != length(lambda))
		{
			stop("Dimension mismatch: Input lambda must contain one number for each network")
		}
		lengths <- unlist(lapply(edge_lengths, sum))
		poisson_parameters <- lambda*lengths
		n.points <- sapply(poisson_parameters, function(x) rpois(n=1, lambda=x))
		return(binomialDesign(n.points)(tree.graphs, edge_lengths, locations, edge_updist, distance_matrices))
	}
	return(replication.function(design.function, replications, rep.variable, rep.values))
}
hardCoreDesign <- function(n, inhibition_region, replications=1, rep.variable = "Time", rep.values)
{
	if(missing(rep.values)) rep.values <- 1:replications
	if(replications != length(rep.values))
	{
		stop("Input rep.values must contain one element for each replication")
	}
	design.function <- function(tree.graphs, edge_lengths, locations, edge_updist, distance_matrices)
	{
		#n can be input as either a vector or a single number
		tmp <- n
		if(length(n) == 1) tmp <- rep(n, length(tree.graphs))
		if(length(inhibition_region) == 1) inhibition_region <- rep(inhibition_region, length(tree.graphs))
		#Check that it has the right number of entries
		if(length(tree.graphs) != length(tmp))
		{
			stop("Dimension mismatch: Input n must contain one number, or one number for each network")
		}
		if(length(tree.graphs) != length(inhibition_region))
		{
			stop("Dimension mismatch: Input inhibition_region must contain one number, or one number for each network")
		}
		initial_points <- binomialDesign(tmp)(tree.graphs, edge_lengths, locations, edge_updist, distance_matrices)
		#Now we are going to go through and strip out points that are closer together than inhibition_region
		final_result <- vector(mode="list", length = length(initial_points))
		cumulative_locID <- 0
		for(i in 1:length(tree.graphs))
		{
			points_this_network <- initial_points[[i]]
			
			#The values returned by uniformDesign have edges indexed uniquely across all the networks, e.g. There is only one edge 1 across all the networks. This is inconvenient, so for now we convert things back so that every network has its own edge number 1. 
			rids <- names(edge_updist[[i]])
			points_this_network$edge <- match(points_this_network$edge, rids)
			
			edges_this_network <- get.edgelist(tree.graphs[[i]])
			edge_lengths_this_network <- edge_lengths[[i]]
			
			distance_matrix_this_network <- distance_matrices[[i]]
			colnames(distance_matrix_this_network) <- rownames(distance_matrix_this_network) <- 1:ncol(distance_matrix_this_network)
						
			point_index <- 1
			distances <- vector(mode="numeric", nrow(points_this_network))
			
			rids_per_point <- points_this_network$edge
			edge_lengths_per_point <- edge_lengths_this_network[rids_per_point]
			
			#distance matrix indexed by rid, or actually by the upstream point of each rid (segment). So we also need the downstream point, so we need to find the line segment one downstream. 
			downstream_rids_per_point <- match(edges_this_network[rids_per_point,1], edges_this_network[, 2])
			
			while(point_index <= nrow(points_this_network))
			{
				current_rid <- points_this_network$edge[point_index]
				downstream_current_rid <- match(edges_this_network[current_rid, 1], edges_this_network[, 2])
				first <- distance_matrix_this_network[downstream_rids_per_point, downstream_current_rid] + edge_lengths_per_point * points_this_network$ratio + edge_lengths_per_point[point_index] * points_this_network$ratio[point_index]
				second <- distance_matrix_this_network[downstream_rids_per_point, current_rid] + edge_lengths_per_point * points_this_network$ratio + edge_lengths_per_point[point_index] * (1 - points_this_network$ratio[point_index])
				third <- distance_matrix_this_network[rids_per_point, downstream_current_rid] + edge_lengths_per_point * (1 - points_this_network$ratio) + edge_lengths_per_point[point_index] * points_this_network$ratio[point_index]
				
				distances <- pmin(first, second, third, na.rm=TRUE)
				subset.bool <- points_this_network$edge == current_rid
				distances[subset.bool] <- edge_lengths_per_point[subset.bool] * abs(points_this_network$ratio[point_index] - points_this_network$ratio[subset.bool])
				
				to.keep <- distances > inhibition_region[i]
				if(any(is.na(to.keep))) stop("Internal error")
				to.keep[point_index] <- TRUE
				
				points_this_network <- points_this_network[to.keep, ]
				rids_per_point <- rids_per_point[to.keep]
				downstream_rids_per_point <- downstream_rids_per_point[to.keep]
				edge_lengths_per_point <- edge_lengths_per_point[to.keep]
				
				point_index <- point_index + 1
			}
			#Here we revert the change mentioned previously, and convert back to having unique edge / rids across all the networks (only one edge number 1, across all the networks). 
			points_this_network$edge <- rids[points_this_network$edge]
			points_this_network$locID <- 1:nrow(points_this_network) + cumulative_locID
			cumulative_locID <- cumulative_locID + nrow(points_this_network)
			
			final_result[[i]] <- points_this_network
		}
		return(final_result)
	}
	return(replication.function(design.function, replications, rep.variable, rep.values))
}
noPoints <- function(tree.graphs, edge_lengths, locations, edge_updist, distance_matrices)
{
	f <- function(x)
	{
		return(data.frame(edge = integer(0), ratio = numeric(0)))
	}
	return(lapply(as.list(1:length(tree.graphs)), f))
}	