get_group_by_individual <- function(association_data, identities = NULL, location = NULL, data_format = c("groups", "individuals")) {

	n_inds <- length(identities)
	if (is.null(identities)) {
		if (data_format == "groups") {
			n_inds <- length(unique(unlist(association_data, use.names = FALSE)))
			ids <- unique(unlist(association_data, use.names = FALSE))
		}
		if (data_format == "individuals") {
			n_inds <- length(unique(association_data[,1]))
			ids <- unique(association_data[,1])
		}
	}
	if (n_inds==0) stop("Error calculating number of individuals")
	
	if (data_format == "groups") {
		if (is.null(names(association_data))) {
			groups <- c(1:length(association_data))
		} else {
			groups <- names(association_data)
		}
	}

	if (data_format == "individuals") {
		groups <- unique(association_data[,2])
		group <- association_data[,2]
	}
	
	if (!is.null(location)) {
		locations <- unique(location)
		if (data_format == "groups") groups <- paste(groups,location,sep="_");
		if (data_format == "individuals") {
			groups <- unique(paste(group,location,sep="_"))
			group <- paste(group,location,sep="_")
		}
	} 

	gbi <- array(0,c(length(groups),n_inds))
	
	colnames(gbi) <- ids
	rownames(gbi) <- groups
	
	if (data_format == "groups") {
		for (i in 1:length(association_data)) { 
			gbi[i,which(ids %in% association_data[[i]])] <- 1
		}
	}
	
	if (data_format == "individuals") {
		for (i in 1:length(groups)) {
			gbi[i,which(ids %in% association_data[which(group %in% groups[i]),1])] <- 1
		}
	}
	
	return(gbi)

}
		