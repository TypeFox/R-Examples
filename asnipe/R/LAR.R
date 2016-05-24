LAR <-
function(group_by_individual, times, timejump, min_time = NULL, max_time = NULL, identities = NULL, which_identities = NULL, locations = NULL, which_locations = NULL, start_time = NULL, end_time = NULL, classes = NULL, which_classes = NULL) {

	
	if (!is.null(colnames(group_by_individual))) { identities <- colnames(group_by_individual) }
	#### CHECK INPUTS
	if (length(dim(group_by_individual)) != 2) { stop("Invalid dimensions for group_by_individual") }
	if (is.null(group_by_individual)) { stop("No group_by_individual data!") }
	if (is.null(times)) { stop("No event time data!") }
	if (is.null(timejump)) { stop("No defined timejump!") }
	if ((length(identities) != ncol(group_by_individual) & !is.null(identities)) == TRUE) { stop("Length of identities does not match number of individuals") }
	if ((length(times) != nrow(group_by_individual) & !is.null(times)) == TRUE) { stop("Length of times does not match number of groups") }
	if ((length(locations) != nrow(group_by_individual) & !is.null(locations)) == TRUE) { stop("Length of locations does not match number of groups") }
	if ((length(classes) != ncol(group_by_individual) & !is.null(classes)) == TRUE) { stop("Length of classes does not match number of individuals") }
	if ((!is.null(which_identities) & is.null(identities)) == TRUE) { stop("Cannot apply which_identities without identities data") }
	if ((!is.null(which_locations) & is.null(locations)) == TRUE) { stop("Cannot apply which_locations without locations data") }
	if ((!is.null(start_time) & is.null(times)) == TRUE) { stop("Cannot apply start_time without times data") }
	if ((!is.null(end_time) & is.null(times)) == TRUE) { stop("Cannot apply end_time without times data") }
	if ((!is.null(which_classes) & is.null(classes)) == TRUE) { stop("Cannot apply which_class without classes data") }

	#### SUBSET THE DATA
	# By identity
	if (!is.null(which_identities)) {
		group_by_individual <- group_by_individual[,which(identities %in% which_identities)]
		if (!is.null(classes)) { classes <- classes[which(identities %in% which_identities)] }
		identities <- identities[which(identities %in% which_identities)]
	}
	
	# By time
	if (!is.null(start_time) & is.null(end_time)) { end_time <- max(times) }
	if (!is.null(end_time) & is.null(start_time)) { start_time <- min(times) }
	if (!is.null(start_time) & !is.null(end_time)) {
		group_by_individual <- group_by_individual[which(times >= start_time & times <= end_time),]
		locations <- locations[which(times >= start_time & times <= end_time)]
		times <- times[which(times >= start_time & times <= end_time)]
	}
	
	# By location
	if (!is.null(which_locations)) {
		group_by_individual <- group_by_individual[which(locations %in% which_locations),]
		locations <- locations[which(locations %in% which_locations)]
		times <- times[which(locations %in% which_locations)]
	}
	
	# By class
	if (!is.null(which_classes)) {
		group_by_individual <- group_by_individual[,which(classes %in% which_classes)]
		if (!is.null(identities)) { identities <- identities[which(classes %in% which_classes)] }
		classes <- classes[which(classes %in% which_classes)]
	}

	#### CALCULATE LAGGED ASSOCIATION RATE
	# Check identities
	if (is.null(identities)) { identities <- c(1:ncol(group_by_individual)) }
	
	# Calculate time steps
	times <- times - min(times)
	if (is.null(min_time)) { min_time <- min(times) }
	if (is.null(max_time)) { max_time <- max(times) }

	timesteps <- seq(min_time,(max_time/timejump),1)
	cat(paste("Timesteps =",length(timesteps)))

	ncells <- ncol(group_by_individual)*(ncol(group_by_individual)-1)
	
	lagged_rates <- matrix(NA,nrow=c(length(timesteps)*ncol(group_by_individual)),ncol=3)
	
	group_by_individual <- group_by_individual[order(times),]
	locations <- locations[order(times)]
	times <- times[order(times)]
	rownames(group_by_individual) <- round(as.numeric(times)/timejump)

	get_associates <- function(i, occur, timesteps, j, identities) {
			reoccur <- occur[as.numeric(rownames(occur)) %in% (as.numeric(rownames(occur))+timesteps[i]),]

			if (!is.null(nrow(occur)) & i!=j) {
			if (nrow(occur)>0) {
				
				if (!is.null(nrow(reoccur))) {
					occur <- occur[which(as.numeric(rownames(occur)) %in% (as.numeric(rownames(reoccur))-timesteps[i])),]
					occur[occur>0] <- 1
					reoccur[reoccur>0] <- 1
					return(cbind(i,sum(reoccur[,-j]*occur[,-j]),sum(reoccur[,j]*occur[,-j])))
				} else {
					occur <- occur[which(rownames(occur)==as.numeric(rownames(occur)[as.numeric(rownames(occur)) %in% (as.numeric(rownames(occur))+timesteps[i])])-timesteps[i]),]
					occur[occur>0] <- 1
					reoccur[reoccur>0] <- 1
					return(cbind(i,sum(reoccur[-j]*occur[-j]),sum(reoccur[j]*occur[-j])))
				}

			}
			} else { return(c(i,NA,NA)) }
	}
	

	for (j in c(1:ncol(group_by_individual))) {
			occur <- group_by_individual[group_by_individual[,j]>0,]
			if (!is.null(nrow(occur))) {
			if (nrow(occur)>0) {
				occur <- t(sapply(by(occur,rownames(occur),colSums),identity))
				if (!is.null(nrow(occur))) {
					lagged_rates[(length(timesteps)*(j-1)+1):(length(timesteps)*j),] <- do.call("rbind",lapply(c(1:length(timesteps)),FUN=get_associates,occur=occur,timesteps=timesteps,j=j, identities=identities))
				}
			} 
			}
	}

	rownames(lagged_rates) <- lagged_rates[,1]
	lagged_rates <- lagged_rates[is.finite(lagged_rates[,3]),]
	lagged_rates <- t(sapply(by(lagged_rates,rownames(lagged_rates),colSums),identity))
	lagged_rates <- lagged_rates[order(as.numeric(rownames(lagged_rates))),c(2,3)]
	return (cbind(log(as.numeric(rownames(lagged_rates))[-1]),lagged_rates[-1,1]/lagged_rates[-1,2]))
}
