# BAMMtools: all functions except for getBAMMlikelihood are undocumented internal functions. 

# This function to compute likelihood exactly as in BAMM.

# break tree into segments, map events. Generate complete matrix for recursive calculations.

# event: 

# index, node, event, timestart, timeend, E0, Et, Dt

# create list of time vectors, of length 1: max number of nodes.
# event_times[[k]] would give the vector of absolute times of events that happened on branch leading to node k
# Event_ID[[k]] would give the index value of the event that happened.


# Generate augmented phylogenetic tree structure with the following components:
#      event_times: a list of length (number of nodes), where event_times[[k]] is 
#                      the vector of absolute times, in order, of events that happened 
#                      on a focal branch. If no event, it is NULL
#      event_id   : a list of length equal to number of nodes, as event_times, but
#                    holding the corresponding event id
#       events    : a dataframe giving parameters and associated nodes (and unique index values)
#                     of the event data.
#      node_event : The event governing the process realized at the node. This will be the first
#                   event encountered as one moves rootwards towards the tips from the focal node
 
eventMatrix <- function(x, phy) {
 
	if (class(x) != 'data.frame') {
		x <- read.csv(x, header = FALSE, stringsAsFactors = FALSE)			
	}

	colnames(x) <- c("generation", "leftchild", "rightchild", "abstime", "lambdainit", "lambdashift", "muinit", "mushift")
	#generation,leftchild,rightchild,abstime,lambdainit,lambdashift,muinit,mushift
	x$index <- 1:nrow(x)
	x$node <- numeric(nrow(x))
	for (i in 1:nrow(x)) {
		if (is.na(x$rightchild[i])) {
			x$node[i] <- which(phy$tip.label == x$leftchild[i])
		} else {
			x$node[i] <- getMRCA(phy, as.character(x[i,2:3]))	
		}
	}
	
	return(x)
}


# A modified version of function eventMatrix(...)
#   but creates a constant-rate birth-death matrix
#   so you can compare likelihoods to diversitree etc.

makeBDeventMatrix <- function(phy, lambda, mu) {
 
	
	xx <- data.frame(generation=1, leftchild=phy$tip.label[1], rightchild = phy$tip.label[length(phy$tip.label)], abstime = 0, lambdainit = lambda, lambdashift = 0, muinit = mu, mushift = 0, stringsAsFactors = FALSE)
	
	colnames(xx) <- c("generation", "leftchild", "rightchild", "abstime", "lambdainit", "lambdashift", "muinit", "mushift")
	#generation,leftchild,rightchild,abstime,lambdainit,lambdashift,muinit,mushift
	xx$index <- 1:nrow(xx)
	xx$node <- numeric(nrow(xx))
	for (i in 1:nrow(xx)){
		if (is.na(xx$rightchild[i])){
			xx$node[i] <- which(phy$tip.label == xx$leftchild[i])
		}else{
			xx$node[i] <- getMRCA(phy, as.character(xx[i,2:3]))	
		}
	}
	
	return(xx)
}

buildTreeWithEventList <- function(phy, events) {
 
	nodeset <- 1:max(phy$edge)
	
	
	phy$event_times <- vector("list", length = length(nodeset))
	phy$event_id <- vector("list", length = length(nodeset))
	phy$node_event <- numeric(length(nodeset))
	
	bt <- branching.times(phy)
	
	phy$events <- events 
	rootnode <- length(phy$tip.label) + 1
	phy <- recursiveAddEventsToTree(phy, rootnode)
	
	
}



recursiveAddEventsToTree <- function(phy, node) {
	
	rootnode <- length(phy$tip.label) + 1
	if (node == rootnode) {
		phy$event_times[[node]] <- 0
		phy$event_id[[node]] <- 1
		phy$node_event[node] <- 1
 		
 		dset <- phy$edge[,2][phy$edge[,1] == node]
		phy <- recursiveAddEventsToTree(phy, dset[1])
		phy <- recursiveAddEventsToTree(phy, dset[2])	
 
	} else {

		parent <- phy$edge[,1][phy$edge[,2] == node]
		events_on_branch <- sum(phy$events$node == node)
		
		if (events_on_branch == 0) {
			phy$node_event[node] <- phy$node_event[parent]
		} else {
			tmp <- phy$events[phy$events$node == node, ]
			tmp <- tmp[order(tmp$abstime), ]
			phy$event_times[[node]] <- tmp$abstime
			phy$event_id[[node]] <- tmp$index
			phy$node_event[node] <- phy$event_id[[node]][events_on_branch]
		}
		
		if (node > length(phy$tip.label)) {
		
			dset <- phy$edge[,2][phy$edge[,1] == node]
			phy <- recursiveAddEventsToTree(phy, dset[1])
			phy <- recursiveAddEventsToTree(phy, dset[2])
		}
 
	}
 	return(phy)
}


# make a list of matrices for each node
# also a vector of E0 and D0 for calculations

# recursively assemble the segment matrix. 
# need seglength value.
# rel_seg is the relative segment length, exactly as specified in BAMM
buildSegmentMatrix <- function(phy, rel_seg){

	nodeset <- 1:max(phy$edge) 			

	phy$D_0 <- rep(1, length(nodeset))
 	phy$E_0 <- numeric(length(nodeset))
	phy$branchsegs <- vector("list", length(nodeset))
 
	bt <- branching.times(phy)
	phy$maxtime <- max(bt)
	
	phy$seg <- rel_seg * phy$maxtime
	
	bt <- max(bt) - bt
	phy$bt <- rep(phy$maxtime, length(nodeset))
	names(phy$bt) <- as.character(nodeset)	
	
	# now includes extant nodes; also speciation times, not branching times
	phy$bt[names(bt)] <- bt

	phy$postorder <- NULL
	phy <- buildSegmentMatricesRecursive(phy, (length(phy$tip.label) + 1))
	
}



buildSegmentMatricesRecursive <- function(phy, node) {

 	seg <- phy$seg
 
 	# these branching times are actually speciation times:
 	#  start at zero, and tip nodes have value equal to max age of tree
	bt <- phy$bt
	
 	if (node > length(phy$tip.label)) {
		dset <- phy$edge[,2][phy$edge[,1] == node]
 		phy <- buildSegmentMatricesRecursive(phy, dset[1])
 		phy <- buildSegmentMatricesRecursive(phy, dset[2])
 	}
	if (node != (length(phy$tip.label) + 1)) {
		#cat("here")
		parent <- phy$edge[,1][phy$edge[,2] == node]
		
		ev <- phy$event_id[[node]]
		et <- phy$event_times[[node]]
		
		branch_end <- bt[as.character(node)]
		branch_start <- bt[as.character(parent)] 
 		start_time <- branch_end
 		end_time <- branch_end		
		
		if (is.null(ev)) {
			
			# If there are no events on branch
			#   do not worry about them
			#   but still do piecewise calculations.
			
			# since ev is null, branch is governed by process of parent node:
 			ev <- phy$node_event[parent]
			
			while (start_time > branch_start) {
				start_time <- start_time - seg
				if (start_time < branch_start) {
					start_time <- branch_start
				}	
				
				tmp <- matrix(c(ev, start_time, end_time, NA, NA, NA, NA, NA), nrow = 1)				
				if (is.null(phy$branchsegs[[node]])){
					phy$branchsegs[[node]] <- tmp
				} else {
					phy$branchsegs[[node]] <- rbind(phy$branchsegs[[node]], tmp)
				}
				
				
				end_time <- start_time
			}
			
 		} else {
			
			while (length(et) > 0) {
				start_time <- start_time - seg
				ll <- length(et)
				if (start_time <= et[ll]){
					start_time <- et[ll]
					et <- et[-ll]				
				}
				tmp <- matrix(c(ev[ll], start_time, end_time, NA, NA, NA, NA, NA), nrow=1)				
				if (is.null(phy$branchsegs[[node]])) {
					phy$branchsegs[[node]] <- tmp
				} else {
					phy$branchsegs[[node]] <- rbind(phy$branchsegs[[node]], tmp)
				}
				
				end_time <- start_time
			}
			
			# now, should be done with all events. Now go through
			#  to beginning of branch
			ev <- phy$node_event[parent]
			
			while (start_time > branch_start) {
				
				
				start_time <- start_time - seg
				if (start_time < branch_start) {
					start_time <- branch_start
				}			 
				tmp <- matrix(c(ev, start_time, end_time, NA, NA, NA, NA, NA), nrow = 1)				
				phy$branchsegs[[node]] <- rbind(phy$branchsegs[[node]], tmp)
				
				end_time <- start_time
				
			}			
			
			
		}
		
 
	}
	
	if (is.null(phy$postorder)) {
		phy$postorder <- node
	} else {
		phy$postorder <- c(phy$postorder, node)
	}
 
 
	return(phy)
}

# likelihood functions for the constant-rate birth-death process
#    used to compute speciation and extinction probabilities
#    on individual branches.

E_func <- function(lam, mu, E0, dt) {
	
	num <- (1 - E0) * (lam - mu);
	denom <- (1 - E0) * lam - (mu - lam * E0) * exp(-(lam - mu) * dt);
	
	return( 1 - num / denom)	
 
}

D_func <- function(lam, mu, E0, D0, dt) {
	
	r <- lam - mu
	num <- (D0 * r^2) * exp((-r) * dt);
	denom <- ( (lam - lam * E0 + exp((-r) * dt) * (lam * E0 - mu)  )  ) ^ 2;
	return(num / denom) 
 
}

############ Likelihood calculation
# Let E_0 and D_0 be initial extinction and data probabilities
# Let E_t and D_t be probabilities computed after some time t, given 
#                initial values.
# Do not track D0 values, always refactor to 1.
# Initialize E0 at tip nodes.
# Loop over postorder sequence
#    for each segment
#        Initialize E_0
#             case (i): if first (tipwards) segment, take E_0 from node. These 
#                       will already have been computed or set.
#             case (ii): otherwise, take previous branch E_t as E_0 
#             case (iii): if finish branch and if identical in state to other branch
#                         set E0 for parent node equal to this value
#             case (iv): if finish branch and not identical in state to other branch
#						 multiply, thus conditioning the extinction probability on the occurrence 
#				         of two lineages at this time
 
#       
#        Compute D_t given this E_0 and add log to likelihood.
#        Compute E_t for segment 
#    At end of segment matrix, E_t for the parent node is set to E_t from 
#                                            the last segment calculation
#         
#       
#
# computeBAMMlikelihood
#     sf:     sampling fraction
#     phy:    phylogenetic tree with all components from buildSegmentMatrixEtc
#
#   Should do this 2 ways for constant-rate process:
#       1. always recomputed E0 for each segment
#       2. do it using segments (the segment advantage is key for time-varying process)

computeBAMMlikelihood <- function(phy, sf = 1, alwaysRecomputeE0 = FALSE, e_prob_condition = "if_different", TOL = 0.001) {
	
	# initial calculation: lets us start with D_0 = 1 at all nodes.
	logLik <- length(phy$tip.label) * log(sf)
	
	phy$E_0[1:length(phy$E_0)] <- -1  # set all values initiall to < 0
	phy$E_0[1:length(phy$tip.label)] <- 1 - sf   # set tip E0
		
	events <- phy$events
	
	 
	for (i in phy$postorder[1:(length(phy$postorder) - 1)]) {
		
		em <- phy$branchsegs[[i]]
		
		for (k in 1:nrow(em)) {
			# elements of branchsegs matrix, in order:
			#    event index, start time , end time, 
			#        E_init for segment, E_final for seg, D_final for seg

			curr_event <- em[k,1]
			index <- which(phy$events$index == curr_event)
			
			# time for start and stop of interval
			#   expressed in units of time since start of current process:
			event_t_start <- em[k,2] -  events$abstime[index]
			event_t_end   <- em[k,3] - events$abstime[index]
			
			lambdainit <- events[index, "lambdainit"]
			lambdashift <- events[index, "lambdashift"]
			muinit <- events[index, "muinit"]
			mushift <- events[index, "mushift"]			
			
			curr_lam <- meanExponentialRate(lambdainit, lambdashift, event_t_start, event_t_end)
			curr_mu <- meanExponentialRate(muinit, mushift, event_t_start, event_t_end)
 	
			tt <- em[k,3] - em[k,2] 
			
			em[k,7] <- curr_lam
			em[k,8] <- curr_mu
			
			if (k == 1 & (!alwaysRecomputeE0)) {
				# if first segment, set to parent node E_init
				em[k,4] <- phy$E_0[i]
			} else {
				em[k,4] <- em[k-1,5]
				# at this point, all E_0 values should be set.	
			}
			
			# compute extinction prob on segment:
			em[k,5] <- E_func(curr_lam, curr_mu, em[k,4], tt)	
			
			# compute speciation prob on segment
			em[k,6] <- D_func(curr_lam, curr_mu, em[k,4], 1.0, tt)
			
			logLik <- logLik + as.numeric( log(em[k,6]) )
	 
		} # for k loop
		if (k == nrow(em)) {
			# at end of branch segments. Ef for last calculation
			# becomes E_0 for the parent node IF identical in state
			# else multiply, thus conditioning on a speciation event.
			
			parent <- phy$edge[,1][phy$edge[,2] == i]
			e0a <- phy$E_0[parent] 
			# value for other descendant branc
			#   if not equal to -1, it will already have been set by other branch
 
			
			if (e0a < 0){
				# value has not been set
				phy$E_0[parent] <- em[nrow(em), 5]
			}else{
				
				# here is value for other lineage
				e0b <-  em[nrow(em), 5] 
 				if (e_prob_condition == "arbitrary"){
					# just making this explicit, taking the "a" branch:
					#cat("arbitrary\n")
					phy$E_0[parent] <- e0a 
				}else if (e_prob_condition == "all_nodes"){
					phy$E_0[parent] <- e0a * e0b
					#cat("all nodes\n")
				}else if (e_prob_condition == "if_different"){
					#cat("if_diff\n")
					delta <- abs(e0a - e0b)
					if (delta > TOL){
						phy$E_0[parent] <- e0a * e0b
					}else{
						phy$E_0[parent] <- e0a 
					}
					
				}else if (e_prob_condition == "random"){
					if (runif(1) < 0.5){
						phy$E_0[parent] <- e0a
					}else{
						phy$E_0[parent] <- e0b
					}
				}else{
					
					stop("Invalid options for e_prob_condition")
				}
			
				# case i: take value abitrarily
				# case ii. condition all nodes
				# case iii. Condition only if extinction probs are different
 				# case iv. take value at random
 				#    e.g., assume 1 lineage is a true parent process
 				
			}
			
			
			#parent <- phy$edge[,1][phy$edge[,2] == i]
			#phy$E_0[parent] <- em[nrow(em), 5]
		}	
		phy$branchsegs[[i]] <- em		
	} # for i loop
 
	# Calculations on nodes:
	# This explicitly conditions on occurrence of a root node
	#    because log(lambda) at the root is not added
	
	nodeset <- phy$postorder[phy$postorder > (length(phy$tip.label) + 1)]
 
	for (i in nodeset) {
		# this block is valid for time-varying rates
		# as exponential speciation rates are computed for each node
		
		curr_lam <- 0
		
		curr <- phy$node_event[i]
		time_from_event <- phy$bt[as.character(i)] - phy$events$abstime[ phy$events$index == curr ]
 
		lam0 <- events[events$index == phy$node_event[i], "lambdainit" ]
		mu <- events[events$index == phy$node_event[i], "muinit" ]
		lshift <- events[events$index == phy$node_event[i], "lambdashift"]
		
		if (lshift <= 0){
			curr_lam <- as.numeric(lam0 * exp(time_from_event * lshift))			
		}else{
			curr_lam <- as.numeric(lam0 * (2 - exp(-lshift * time_from_event)))
		}
		

		logLik <- logLik + log(curr_lam)
	
	}
 
 
 	## To condition on survival:
 	#  Need the extinction probs of process at 2 basal branches:
	
	dset <- phy$edge[phy$edge[,1] == (length(phy$tip.label) + 1), 2]
	e1 <- phy$branchsegs[[dset[1]]]
	e2 <- phy$branchsegs[[dset[2]]]
 
 	
 	E_prob_1  <- e1[nrow(e1), 5]
 	E_prob_2 <- e2[nrow(e2), 5]
 	# These are the probability that a single lineage
 	#    does not go extinct.
 	#  So probability that both basal branches persist 
 	#  is (1 - p(extinct))^2
 	#  This is exactly isometric with the diversitree likelihood
  	#  
  	# But treat R and L branches separately as may have different values
  
 	logLik <- logLik - log(1 - E_prob_1) - log(1 - E_prob_2)
	
	phy$logLik <- logLik
	
	return(phy)
	
}
 
meanExponentialRate <- function(rate_init, rate_shift, t_start, t_end) {
	
	delta_T = t_end - t_start;
    integrated = 0.0;
    
    if (rate_shift < 0) {
        integrated = (rate_init / rate_shift) * (exp(rate_shift * t_end) - exp(rate_shift * t_start));
    } else if (rate_shift > 0) {
        integrated = rate_init * (2 * delta_T + (1.0 / rate_shift) *
                       (exp(-rate_shift * t_end) - exp(-rate_shift * t_start)));
    } else {
        integrated = rate_init * delta_T;
    }
	
	return(integrated / delta_T)
} 


BAMMlikelihood <- function(phy, eventdata, gen = 'last', segLength = 0.02, sf = 1, return.intermediates = FALSE, e_prob_condition = "if_different", ...) {
	#gen can be a number 1 -> numberOfGenerations, or 'last', or 'all'
	#segLength is segLength value used in BAMM
	#sf is sampling fraction
	
	if (sum(names(eventdata) %in% c('lambda', 'mu')) == 2 & length(eventdata == 2)){
		eventdata <- makeBDeventMatrix(phy, eventdata['lambda'], eventdata['mu'])
		
	}
	
	if (class(phy) == 'character') {
		phy <- read.tree(phy)
	}
	
	if (class(eventdata) != 'data.frame') {
		eventdatafile <- eventdata
		eventdata <- read.csv(eventdata, header = FALSE, stringsAsFactors = FALSE)			
	}
	
	colnames(eventdata) <- c("generation", "leftchild", "rightchild", "abstime", "lambdainit", "lambdashift", "muinit", "mushift")
	#if eventdata already had header, fix
	if (eventdata[1,1] == 'generation') {
		eventdata <- read.csv(eventdatafile, header = TRUE, stringsAsFactors = FALSE)
	}
	
	#check that supplied generation number is valid
	if (is.numeric(gen)) {
		if (any(!gen %in% eventdata$generation)) {
			stop('Supplied generation number is not valid.')
		}
	}
	
	if (gen == 'last') {
		gen <- tail(eventdata$generation, 1)
	}
	
	if (gen == 'all') {
		gen <- unique(eventdata$generation)
	}
	
	#extract requested generation
	if (length(gen) != length(unique(eventdata$generation))) {
		eventdata <- eventdata[sort(unlist(sapply(gen, function(x) which(eventdata$generation == x)))),]
	}
	
	eventList <- split(eventdata, eventdata$generation)
	
	if (return.intermediates){
		res <- vector(length=length(eventList), mode = "list")
		for (i in 1:length(eventList)) {
		 
			eMat <- eventMatrix(eventList[[i]], phy)
			phy2 <- buildTreeWithEventList(phy, eMat)
			phy2 <- buildSegmentMatrix(phy2, segLength)	
			res[[i]] <- computeBAMMlikelihood(phy2, sf, e_prob_condition = e_prob_condition, ...)
		}		
		return(res)		
	}else{
		res <- vector(length = length(eventList))
		for (i in 1:length(eventList)) {
		 
			eMat <- eventMatrix(eventList[[i]], phy)
			phy2 <- buildTreeWithEventList(phy, eMat)
			phy2 <- buildSegmentMatrix(phy2, segLength)	
			res[i] <- computeBAMMlikelihood(phy2, sf, e_prob_condition = e_prob_condition, ...)$logLik
		}		
		return(res)
	}
 
}


 

# x is event data frame with locations 
# parameters (lambda, mu etc) are not used.
optimize1shiftBAMM <- function( phy, x , rel_seg = 1, e_prob_condition = "if_different"){
	
	emat <- eventMatrix(x, phy)
	phy2 <- buildTreeWithEventList(phy, emat)
	phy2 <- buildSegmentMatrix(phy2, rel_seg)
	
	lfx <- function(x){
		phy2$events$lambdainit[1:2] <- exp(x[1:2])
		phy2$events$muinit[1:2] <- exp(x[3:4])
		phy2$events$lambdashift[1:2] <- c(0,0)
		phy2$events$mushift[1:2] <- c(0,0)

		tmp <- computeBAMMlikelihood(phy2, sf=1, e_prob_condition = e_prob_condition) 
		return(tmp$logLik)
	}
	
	# usually right here I would sample random parameters:
	# like runif(2, 0, 0.5) for lambda
	#  then eps <- runif(2, 0, 1) for relative extinction
	#  then back-compute the extinction rates by multiplying eps * lambda
	
	# the danger is that you need a try statement to catch 
	# if the function cannot be evaluated at the initial parameters
	
	#pars_init <- c(0.1, 0.1, 0.01, 0.01)
	
	lam_init <- runif(2, 0, 0.5)
	mu_init <- runif(2, 0, 1) * lam_init
	pars_init <- c(lam_init, mu_init)
	
	res <- optim(log(pars_init), fn=lfx, method = "Nelder", control=list(fnscale = -1))
 
	rr <- list(logLik = res$value, lambda = exp(res$par[1:2]), mu = exp(res$par[3:4]), conv=res$convergence, e_prob_condition = e_prob_condition)
	return(rr)
}











