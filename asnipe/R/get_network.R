get_network <- function(association_data, data_format = "GBI", association_index = "SRI", identities = NULL, which_identities = NULL, times = NULL, locations = NULL, which_locations = NULL, start_time = NULL, end_time = NULL, classes = NULL, which_classes = NULL) {

	#### CHECK INPUTS
	if (is.null(association_data)) { stop("No association_data data!") }
	if (length(dim(association_data)) != 2 & data_format=="GBI") { stop("Invalid dimensions for association_data") }
	if (length(dim(association_data)) != 3 & data_format=="SP") { stop("Invalid dimensions for association_data") }
	if ((length(identities) != ncol(association_data) & !is.null(identities)) == TRUE) { stop("Length of identities does not match number of individuals") }
	if ((length(times) != nrow(association_data) & !is.null(times)) == TRUE) { stop("Length of times does not match number of groups") }
	if ((length(locations) != nrow(association_data) & !is.null(locations)) == TRUE) { stop("Length of locations does not match number of groups") }
	if ((length(classes) != ncol(association_data) & !is.null(classes)) == TRUE) { stop("Length of classes does not match number of individuals") }
	if ((!is.null(which_identities) & is.null(identities)) == TRUE) { stop("Cannot apply which_identities without identities data") }
	if ((!is.null(which_locations) & is.null(locations)) == TRUE) { stop("Cannot apply which_locations without locations data") }
	if ((!is.null(start_time) & is.null(times)) == TRUE) { stop("Cannot apply start_time without times data") }
	if ((!is.null(end_time) & is.null(times)) == TRUE) { stop("Cannot apply end_time without times data") }
	if ((!is.null(which_classes) & is.null(classes)) == TRUE) { stop("Cannot apply which_class without classes data") }
	if (!any(association_index %in% c("SRI","HWI"))) { stop("Unknown association_index") }

	

	#### SUBSET THE DATA
	# By identity
	if (!is.null(which_identities)) {
		if (data_format=="GBI") association_data <- association_data[,which(identities %in% which_identities)]
		if (data_format=="SP") association_data <- association_data[,which(identities %in% which_identities),which(identities %in% which_identities)]
		identities <- identities[which(identities %in% which_identities)]
	}
	
	# By time
	if (!is.null(start_time) & is.null(end_time)) { end_time <- max(times) }
	if (!is.null(end_time) & is.null(start_time)) { start_time <- min(times) }
	if (!is.null(start_time) & !is.null(end_time)) {
		subs <- which(times >= start_time & times <= end_time)
		if (data_format=="GBI") association_data <- association_data[subs,]
		if (data_format=="SP") association_data <- association_data[subs,,]
		locations <- locations[subs]
		times <- times[subs]
	}
	
	# By location
	if (!is.null(which_locations)) {
		subs <- which(locations %in% which_locations)
		if (data_format=="GBI") association_data <- association_data[subs,]
		if (data_format=="SP") association_data <- association_data[subs,,]
		locations <- locations[subs]
		times <- times[subs]
	}
	
	# By class
	if (!is.null(which_classes)) {
		if (data_format=="GBI") association_data <-  association_data[,which(classes %in% which_classes)]
		if (data_format=="SP") association_data <- association_data[,which(classes %in% which_classes),which(classes %in% which_classes)]
		identities <- identities[which(classes %in% which_classes)]
	}
	
	#### GENERATE NETWORK
	### Calculate Network
	do.SR <- function(GroupBy,input,association_index){
		jumps <- c(seq(0,ncol(input),50))
		if (max(jumps) < ncol(input)) {
			jumps <- c(jumps,ncol(input))
		}
		out <- matrix(nrow=0,ncol=1)
		for (i in 1:(length(jumps)-1)) {
			tmp <- input[ ,GroupBy] + input[,(jumps[i]+1):jumps[i+1]]
			if (length(tmp) > nrow(input)) {
				x <- colSums(tmp==2)
			} else {
				x <- sum(tmp==2)
			}
			if (length(tmp) > nrow(input)) {
				yab <- colSums(tmp==1)
			} else {
				yab <- sum(tmp==1)
			}
			if (association_index == "SRI") {
				out <- c(out, x / (x + yab))
			} else if (association_index == "HWI") {
				out <- c(out, x / (x + 0.5*yab))
			}
		}
		out
	}
	
	do.SR2 <- function (i, a,association_index) {
		# how many times 1 seen together with all others
		x <- apply(a[,i,],2,sum)

		# how many times 1 but not others in a sampling period and vice versa
		n <- apply(a,1,rowSums)
		n[n>0] <- 1
		seen <- t(apply(n,1,function(x) x-n[i,]))
		ya <- rowSums(seen<0)
		yb <- rowSums(seen>0)

		# how many times 1 and others seen but not together
		seen <- t(apply(n,1,function(x) x+n[i,]))
		yab <- rowSums(seen>1) - x
		
		if (association_index == "SRI") {
			out <- x / (x + ya + yb + yab)
		} else if (association_index == "HWI") {
			out <- x / (x + ya + yb + 0.5*yab)
		}
		return(out)
	}
	
	cat(paste("Generating ", ncol(association_data), " x ", ncol(association_data), " matrix\n"))

	if (data_format=="GBI") fradj_sorted <- do.call("rbind",lapply(seq(1,ncol(association_data),1),FUN=do.SR,input=association_data, association_index))
	
	if (data_format=="SP") fradj_sorted <- do.call("rbind",lapply(seq(1,ncol(association_data),1),FUN=do.SR2,a=association_data, association_index))
		
	
	fradj_sorted[is.nan(fradj_sorted)] <- 0
	diag(fradj_sorted) <- 0

	if (!is.null(identities)) {
		colnames(fradj_sorted) <- identities
		rownames(fradj_sorted) <- identities
	} else if (!is.null(colnames(association_data))) {
		colnames(fradj_sorted) <- colnames(association_data)
		rownames(fradj_sorted) <- colnames(association_data)
	}
	
	return(fradj_sorted)
}
