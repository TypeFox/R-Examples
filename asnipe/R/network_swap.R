network_swap <- function(association_data, data_format = "GBI", swaps=1, association_index = "SRI", association_matrix = NULL, locations = NULL, classes = NULL, days = NULL, within_day = FALSE, within_location = FALSE, within_class = FALSE) {
	
	permutations <- swaps
	
	#### CHECK INPUTS
	if (is.null(association_data)) { stop("No association_data data!") }
	if (length(dim(association_data)) != 2 & data_format=="GBI") { stop("Invalid dimensions for association_data") }
	if (length(dim(association_data)) != 3 & data_format=="SP") { stop("Invalid dimensions for association_data") }
	if ((length(locations) != nrow(association_data) & !is.null(locations)) == TRUE) { stop("Length of locations does not match number of groups") }
	if ((length(classes) != ncol(association_data) & !is.null(classes)) == TRUE) { stop("Length of classes does not match number of individuals") }
	if (!any(association_index %in% c("SRI","HWI"))) { stop("Unknown association_index") }
	if (within_day & is.null(days)) { stop("Cannot constrict within days if days are not supplied") }
	if (within_location & is.null(locations)) { stop("Cannot constrict within location if locations are not supplied") }
	if (within_class & is.null(classes)) { stop("Cannot constrict within class if classes are not supplied") }
	
	if (!within_day) { days <- rep(1,nrow(association_data)) }
	if (!within_location) { locations <- rep(1,nrow(association_data)) }
	if (!within_class) { classes <- rep(1,ncol(association_data)) }
	
	
	#### GENERATE NETWORK IF REQUIRED
	### Calculate Network
		
	if (is.null(association_matrix)) {
	
		fradj_sorted <- get_network(association_data,data_format)

	} else {
		fradj_sorted <- association_matrix
	}


	#### DO PERMUTATIONS
	association_data_perm <- association_data
	fradj_sorted2 <- fradj_sorted
	n_inds <- ncol(association_data_perm)
	
	# Calculate network
	do.SR_perm <- function(GroupBy,input){
			tmp <- input[ ,GroupBy] + input
			x <- colSums(tmp==2)
			yab <- colSums(tmp==1)
			out <- (x / (x + yab))
		out
	}
	
	do.SR2_perm <- function(i, a, association_index) {
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

	count <- 1
	##  GET PERMUTATION MATRICES
	for (n in c(1:permutations)) {
		repeat {
			a <- which(association_data_perm>0,arr.ind=TRUE)
			s <- sample(1:nrow(a),1)
			first <- a[s,]
			if (data_format=="GBI") {
				a <- a[which(locations[as.numeric(a[,1])] == locations[as.numeric(a[s,1])] & days[as.numeric(a[,1])] == days[as.numeric(a[s,1])] & classes[as.numeric(a[,2])] == classes[as.numeric(a[s,2])]),]
				second <- a[sample(1:nrow(a),1),]
				if (first[1]!=second[1]&first[2]!=second[2]&(sum(association_data_perm[first[1],first[2]]) > 0 & sum(association_data_perm[second[1],second[2]]) > 0) & (sum(association_data_perm[second[1],first[2]]) == 0 & sum(association_data_perm[first[1],second[2]]) == 0)) { break; }
			}
			if (data_format=="SP") {
				a <- a[which(a[,1] == a[s,1] & classes[as.numeric(a[,2])] == classes[as.numeric(a[s,2])] & classes[as.numeric(a[,3])] == classes[as.numeric(a[s,3])]),]
				second <- a[sample(1:nrow(a),1),]
				if (first[2]!=second[2]&first[3]!=second[3]&(sum(association_data_perm[first[1],first[2],first[3]]) > 0 & sum(association_data_perm[second[1],second[2],second[3]]) > 0) & (sum(association_data_perm[second[1],second[2],first[3]]) == 0 & sum(association_data_perm[first[1],first[2],second[3]]) == 0)) { break; }
			}
				
		}
		
		if (data_format=="GBI") {
			association_data_perm[second[1],first[2]] <- association_data_perm[first[1],first[2]]
			association_data_perm[first[1],second[2]] <- association_data_perm[second[1],second[2]]
			association_data_perm[first[1],first[2]] <- 0
			association_data_perm[second[1],second[2]] <- 0
		
			tmp1 <- do.SR_perm(first[2],association_data_perm)
			tmp2 <- do.SR_perm(second[2],association_data_perm)
		}
		if (data_format=="SP") {
			association_data_perm[first[1],second[2],first[3]] <- association_data_perm[first[1],first[2],first[3]]
			association_data_perm[first[1],first[2],second[3]] <- association_data_perm[first[1],second[2],second[3]]
			association_data_perm[first[1],first[2],first[3]] <- 0
			association_data_perm[first[1],second[2],second[3]] <- 0
		
			tmp1 <- do.SR2_perm(first[2],association_data_perm, association_index)
			tmp2 <- do.SR2_perm(second[2],association_data_perm, association_index)
		}

		fradj_sorted2[,first[2]] <- tmp1
		fradj_sorted2[first[2],] <- tmp1
		fradj_sorted2[,second[2]] <- tmp2
		fradj_sorted2[second[2],] <- tmp2
		
		diag(fradj_sorted2) <- 0
		
	}
	
	return(list(Association_index=fradj_sorted2, Association_data=association_data_perm))
}
