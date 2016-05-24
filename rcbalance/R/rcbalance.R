rcbalance <-
function(distance.structure, near.exact = NULL, fb.list = NULL, treated.info = NULL, control.info = NULL, exclude.treated = FALSE, target.group = NULL,  k = 1, penalty = 3, tol = 1e-5){
		
####################  CHECK INPUT #################### 
	if(!(k > 0)){
		stop("k must be larger than zero")
	}
	if(!(penalty > 1)){
		stop("penalty argument must be greater than 1")
	}
	#exclude.treated is incompatible with k > 1 and with target distributions other than default treated distribution
	if(exclude.treated && (!is.null(target.group) || k > 1)){
		stop("exclude.treated = TRUE is not compatible with target.group arguments or k arguments greater than 1")
	}

	if(!is.null(near.exact)){
		if(!(!is.null(treated.info) && !is.null(control.info))){
			stop('treated.info and control.info must both be specified when near.exact is used')
		}
		if(ncol(treated.info) != ncol(control.info)){
			stop('treated.info and control.info must have the same number of columns')
		}
		if(!all(colnames(treated.info) == colnames(control.info))){
			stop('treated.info and control.info must have identical column names')
		}	
		treated.control.info <- rbind(treated.info, control.info)
		if(class(distance.structure) %in% c('matrix', 'InfinitySparseMatrix', 'BlockedInfinitySparseMatrix')){		
			if(nrow(treated.info) != nrow(distance.structure)){
				stop('Dimensions of treated.info and distance.structure do not agree')
			}
			if(nrow(control.info) != ncol(distance.structure)){
				stop('Dimensions of control.info and distance.structure do not agree')				
			}
		}else{
			#make sure number of treated in distance.structure and treated.info agree
			if(nrow(treated.info) != length(distance.structure)){
				stop('treated.info and distance.structure specify different numbers of treated units')
			}
			#make sure number of controls in distance.structure and control.info agree, i.e. max control index in each list element must not exceed row count in matrix of all subjects
			if(nrow(control.info) < max(laply(distance.structure, function(x) max(c(as.numeric(names(x)),0))))){
				stop('Not all control units in distance.structure have information in control.info')
			}		
		}	
		if(!(all(near.exact %in% colnames(treated.info)))){
			stop('near.exact contains variable names not present in colnames(treated.info)')
		}
	}
	
	if(!is.null(fb.list)){
		if(is.null(target.group) && !is.null(near.exact)){
			#don't need to repeat checks 
 			target.group <- treated.info
			target.control.info <- treated.control.info
		 }else{ 
		 	if(is.null(target.group)){
		 		target.group <- treated.info
		 	}
			if(!(!is.null(treated.info) && !is.null(control.info))){
				stop('treated.info and control.info must both be specified when fb.list is used')
			}
			if(ncol(treated.info) != ncol(control.info)){
				stop('treated.info and control.info must have the same number of columns')
			}
			if(!all(colnames(treated.info) == colnames(control.info))){
				stop('treated.info and control.info must have identical column names')
			}	
			target.control.info <- rbind(target.group, control.info)
			if(class(distance.structure) %in% c('matrix', 'InfinitySparseMatrix', 'BlockedInfinitySparseMatrix')){		
				if(nrow(target.group) != nrow(distance.structure)){
					stop('target.group and distance.structure dimensions do not agree')
				}
				if(nrow(control.info) != ncol(distance.structure)){
					stop('control.info and distance.structure dimensions do not agree')					
				}
			}else{
				#make sure number of treated in distance.structure and target.group agree
				if(nrow(target.group) != length(distance.structure)){
				stop('target.group and distance.structure specify different numbers of treated units')
				}
				#make sure number of controls in distance.structure and control.info agree, i.e. max control index in each list element must not exceed row count in matrix of all subjects
				if(nrow(control.info) < max(laply(distance.structure, function(x) max(c(as.numeric(names(x)),0))))){
					stop('Not all control units in distance.structure have information in control.info')
				}		
			}
			if(!all(unlist(fb.list) %in% colnames(target.group))){
				stop('fb.list contains variable names not given in column names of other arguments')
			} 
			if(length(fb.list) > 1){
				for(i in c(1:(length(fb.list)-1)))	if(!all(fb.list[[i]] %in% fb.list[[i+1]])){
					stop('Elements of fb.list must contain all variables listed in previous elements')
				}	
			}				
		}	
	}


######## SET UP TREATED-CONTROL PORTION OF NETWORK	#########
	if (inherits(distance.structure, c('matrix', 'InfinitySparseMatrix'))) {
		match.network <- dist2net.matrix(distance.structure,k, exclude.treated = exclude.treated)
	} else if (!is.null(fb.list)){ #specify number of controls
		match.network <- dist2net(distance.structure,k, exclude.treated = exclude.treated, ncontrol = nrow(control.info))		
	} else {
		match.network <- dist2net(distance.structure,k, exclude.treated = exclude.treated)
	}
	
####################  ADD FINE BALANCE CONSTRAINTS #################### 

	if(!is.null(fb.list)){
		for(my.layer in fb.list){
				interact.factor <- apply(target.control.info[,match(my.layer, colnames(target.control.info)), drop = FALSE],1, function(x) paste(x, collapse ='.'))
				match.network <- add.layer(match.network, interact.factor)
		}		
	}
		
	match.network <- penalty.update(match.network, newtheta = penalty) 

#################### ADD NEAR EXACT PENALTIES ######################
	
	if(!is.null(near.exact)){
		interact.factor <- apply(treated.control.info[,match(near.exact, colnames(treated.control.info)), drop = FALSE],1, function(x) paste(x, collapse ='.'))
		match.network <- penalize.near.exact(match.network, interact.factor)	
	}


############################ RUN MATCH ##############################
    #convert costs to integers if necessary
    #h/t optmatch developers for ideas about how to do this nicely 
	cost <- match.network$cost
    if(any(cost != round(cost))){
    	intcost <- round(cost/tol)
    	#use smaller scaling factor if possible
    	searchtol <- 10^(-c(1:floor(log10(.Machine$integer.max))))
    	searchtol <- searchtol[searchtol > tol]
    	for (newtol in searchtol){
    		new.intcost <- round(intcost*tol/newtol)
    		if (any(new.intcost != intcost*tol/newtol)) break
    		tol <- newtol
    		intcost <- new.intcost
    	}
    	cost <- intcost    		
	}
	if (any(is.na(as.integer(cost)))) {
		stop('Integer overflow in penalties!  Run with a higher tolerance, a lower penalty value, or fewer levels of fine balance.')
	} 	
	o <- callrelax(match.network)	
	if(o$feasible == 0){
		stub.message <- 'Match is infeasible or penalties are too large for RELAX to process! Consider reducing penalty or raising tolerance'
		if(k > 1){	
			#print()
			stop(paste(stub.message, 'or reducing k.'))
		}
		if(!exclude.treated){
			#print()
			stop(paste(stub.message, 'or setting exclude.treated = TRUE.'))
		}
		#print()
		stop(paste(stub.message, '.', sep =''))
	}
	
	
	#################### PREPARE OUTPUT #################### 	
	#make a |T| x k matrix with rownames equal to index of treated unit and indices of its matched controls stored in each row
	x <- o$x[1:match.network$tcarcs]	
	match.df <- data.frame('treat' = as.factor(match.network$startn[1:match.network$tcarcs]), 'x' = x, 'control' = match.network$endn[1:match.network$tcarcs])
	matched.or.not <- daply(match.df, .(match.df$treat), function(treat.edges) c(as.numeric(as.character(treat.edges$treat[1])), sum(treat.edges$x)), .drop_o = FALSE)
	if(any(matched.or.not[,2] == 0)){
		match.df <- match.df[-which(match.df$treat %in% matched.or.not[which(matched.or.not[,2] == 0),1]),]
	}
	match.df$treat <- as.factor(as.character(match.df$treat))
	matches <- as.matrix(daply(match.df, .(match.df$treat), function(treat.edges) treat.edges$control[treat.edges$x == 1], .drop_o = FALSE))

	#make a contingency table for each fine balance factor 
	if(is.null(fb.list)){
		fb.tables <- NULL
	}else{
		#variables for matched subjects only
		matched.info <- rbind(treated.info[as.numeric(rownames(matches)),,drop = FALSE], control.info[as.vector(matches) - sum(match.network$z),,drop = FALSE])
		treatment.status <- c(rep(1, nrow(matches)), rep(0, k*nrow(matches)))
		#for each fine balance level k, make a vector of nu_k values for the matched subjects	
		interact.factors.matched = llply(fb.list, function(my.layer) as.factor(apply(matched.info[,match(my.layer, colnames(matched.info)), drop = FALSE],1, function(x) paste(x, collapse ='.'))))
		fb.tables <- llply(interact.factors.matched, function(inter.fact) table('balance.variable' = inter.fact, treatment.status))	
	}
	
	#need to decrement match indices to ensure controls are numbered 1:nc again
	final.matches <- matrix(matches - sum(match.network$z), ncol =k, dimnames = list(rownames(matches),1:k))
	final.matches <- final.matches[order(as.numeric(rownames(final.matches))), , drop=FALSE]
	return(list('matches' = final.matches, 'fb.tables' = fb.tables))
}
