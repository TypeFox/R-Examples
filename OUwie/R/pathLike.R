#RECONSTRUCTION OF EVOLUTIONARY PATHWAYS BASED ON FARRIS (1973) and DESCRIBED BY STEEL AND PENNY (2000):

#written by Jeremy M. Beaulieu

PathLik <- function(phy, data, p, hrm=FALSE, rate.cat, ntraits=NULL, charnum=1, rate.mat=NULL, model=c("ER", "SYM", "ARD"), root.p=NULL, include.maps=TRUE){
	#Note: Does not like zero branches at the tips. Here I extend these branches by just a bit:
	phy$edge.length[phy$edge.length<=1e-5]=1e-5
	if(hrm==FALSE){
		if(ntraits==1){
			data.sort<-data.frame(data[,charnum+1],data[,charnum+1],row.names=data[,1])
		}
		if(ntraits==2){
			data.sort<-data.frame(data[,2], data[,3], row.names=data[,1])
		}
		if(ntraits==3){
			data.sort<-data.frame(data[,2], data[,3], data[,4], row.names=data[,1])
		}
	}
	else{
		data.sort <- data.frame(data[,2], data[,2],row.names=data[,1])
	}
	data.sort<-data.sort[phy$tip.label,]
	#Some initial values for use later
	k=2
	obj <- NULL
	nb.tip <- length(phy$tip.label)
	nb.node <- phy$Nnode
	#Builds the rate matrix based on the specified rate.cat. Not exactly the best way
	#to go about this, but it is the best I can do for now -- it works, so what me worry?
	if(hrm==TRUE){
		if(is.null(rate.mat)){
			rate <- rate.mat.maker(hrm=TRUE,rate.cat=rate.cat)
			rate[is.na(rate)] <- max(rate,na.rm=TRUE)+1
		}
		else{
			rate <- rate.mat
			rate[is.na(rate)] <- max(rate,na.rm=TRUE)+1
		}
		#Makes a matrix of tip states and empty cells corresponding 
		#to ancestral nodes during the optimization process.	
		x <- data.sort[,1]
		TIPS <- 1:nb.tip
		
		for(i in 1:nb.tip){
			if(is.na(x[i])){x[i]=2}
		}
		if (rate.cat == 1){
			liks <- matrix(0, nb.tip + nb.node, k*rate.cat)
			TIPS <- 1:nb.tip
			for(i in 1:nb.tip){
				if(x[i]==0){liks[i,1]=1}
				if(x[i]==1){liks[i,2]=1}
				if(x[i]==2){liks[i,1:2]=1}
			}
		}
		if (rate.cat == 2){
			liks <- matrix(0, nb.tip + nb.node, k*rate.cat)
			for(i in 1:nb.tip){
				if(x[i]==0){liks[i,c(1,3)]=1}
				if(x[i]==1){liks[i,c(2,4)]=1}
				if(x[i]==2){liks[i,1:4]=1}
			}
		}
		if (rate.cat == 3){
			liks <- matrix(0, nb.tip + nb.node, k*rate.cat)
			for(i in 1:nb.tip){
				if(x[i]==0){liks[i,c(1,3,5)]=1}
				if(x[i]==1){liks[i,c(2,4,6)]=1}
				if(x[i]==2){liks[i,1:6]=1}
			}
		}
		if (rate.cat == 4){
			liks <- matrix(0, nb.tip + nb.node, k*rate.cat)
			for(i in 1:nb.tip){
				if(x[i]==0){liks[i,c(1,3,5,7)]=1}
				if(x[i]==1){liks[i,c(2,4,6,8)]=1}
				if(x[i]==2){liks[i,1:8]=1}
			}
		}
		if (rate.cat == 5){
			liks <- matrix(0, nb.tip + nb.node, k*rate.cat)
			for(i in 1:nb.tip){
				if(x[i]==0){liks[i,c(1,3,5,7,9)]=1}
				if(x[i]==1){liks[i,c(2,4,6,8,10)]=1}
				if(x[i]==2){liks[i,1:10]=1}
			}
		}
		if (rate.cat == 6){
			liks <- matrix(0, nb.tip + nb.node, k*rate.cat)
			for(i in 1:nb.tip){
				if(x[i]==0){liks[i,c(1,3,5,7,9,11)]=1}
				if(x[i]==1){liks[i,c(2,4,6,8,10,12)]=1}
				if(x[i]==2){liks[i,1:12]=1}
			}
		}
		if (rate.cat == 7){
			liks <- matrix(0, nb.tip + nb.node, k*rate.cat)
			for(i in 1:nb.tip){
				if(x[i]==0){liks[i,c(1,3,5,7,9,11,13)]=1}
				if(x[i]==1){liks[i,c(2,4,6,8,10,12,14)]=1}
				if(x[i]==2){liks[i,1:14]=1}
			}
		}
		Q <- matrix(0, k*rate.cat, k*rate.cat)
		tranQ <- matrix(0, k*rate.cat, k*rate.cat)
	}
	if(hrm==FALSE){
		#Imported from Jeffs rayDISC -- will clean up later, but for now, it works fine:
		if(ntraits == 1){
			k <- 1
			factored <- factorData(data.sort) # was acting on data, not data.sort			
			nl <- ncol(factored)
			obj <- NULL
			nb.tip<-length(phy$tip.label)
			nb.node <- phy$Nnode
			
			if(is.null(rate.mat)){
				rate <- rate.mat.maker(hrm=FALSE,ntraits=ntraits,nstates=nl,model=model)
				rate[is.na(rate)] <- max(rate,na.rm=TRUE)+1
			}
			else{
				rate<-rate.mat
				rate[is.na(rate)] <- max(rate,na.rm=TRUE)+1
			}
			
			stateTable <- NULL # will hold 0s and 1s for likelihoods of each state at tip
			for(column in 1:nl){
				stateTable <- cbind(stateTable,factored[,column])
			}
			colnames(stateTable) <- colnames(factored)
			
			ancestral <- matrix(0,nb.node,nl) # all likelihoods at ancestral nodes will be 0
			liks <- rbind(stateTable,ancestral) # combine tip likelihoods & ancestral likelihoods
			rownames(liks) <- NULL
		}
		if(ntraits==2){
			k=2
			nl=2
			if(is.null(rate.mat)){
				rate <- rate.mat.maker(hrm=FALSE,ntraits=ntraits,model=model)
				rate[is.na(rate)] <- max(rate,na.rm=TRUE)+1
			}
			else{
				rate <- rate.mat
				rate[is.na(rate)] <- max(rate,na.rm=TRUE)+1
			}
			x<-data.sort[,1]
			y<-data.sort[,2]
			
			liks <- matrix(0, nb.tip + nb.node, nl^k)
			TIPS <- 1:nb.tip
			for(i in 1:nb.tip){
				if(is.na(x[i])){
					x[i]=2
				}
				if(is.na(y[i])){
					y[i]=2
				}
			}
			for(i in 1:nb.tip){
				if(x[i]==0 & y[i]==0){liks[i,1]=1}
				if(x[i]==0 & y[i]==1){liks[i,2]=1}
				if(x[i]==1 & y[i]==0){liks[i,3]=1}
				if(x[i]==1 & y[i]==1){liks[i,4]=1}
				if(x[i]==2 & y[i]==0){liks[i,c(1,3)]=1}
				if(x[i]==2 & y[i]==1){liks[i,c(2,4)]=1}
				if(x[i]==0 & y[i]==2){liks[i,c(1,2)]=1}
				if(x[i]==1 & y[i]==2){liks[i,c(3,4)]=1}
				if(x[i]==2 & y[i]==2){liks[i,1:4]=1}
			}
		}
		if(ntraits==3){
			k = 3
			nl = 2
			if(is.null(rate.mat)){
				rate <- rate.mat.maker(hrm=FALSE,ntraits=ntraits,model=model)
				rate[is.na(rate)]<-max(rate,na.rm=TRUE)+1
			}
			else{
				rate <- rate.mat
				rate[is.na(rate)] <- max(rate,na.rm=TRUE)+1
			}			
			x <- data.sort[,1]
			y <- data.sort[,2]
			z <- data.sort[,3]
			
			liks <- matrix(0, nb.tip + nb.node, nl^k)
			TIPS <- 1:nb.tip
			for(i in 1:nb.tip){
				if(is.na(x[i])){
					x[i]=2
				}
				if(is.na(y[i])){
					y[i]=2
				}
				if(is.na(z[i])){
					z[i]=2
				}
			}
			for(i in 1:nb.tip){
				if(x[i]==0 & y[i]==0 & z[i]==0){liks[i,1]=1}
				if(x[i]==1 & y[i]==0 & z[i]==0){liks[i,2]=1}
				if(x[i]==0 & y[i]==1 & z[i]==0){liks[i,3]=1}
				if(x[i]==0 & y[i]==0 & z[i]==1){liks[i,4]=1}
				if(x[i]==1 & y[i]==1 & z[i]==0){liks[i,5]=1}
				if(x[i]==1 & y[i]==0 & z[i]==1){liks[i,6]=1}
				if(x[i]==0 & y[i]==1 & z[i]==1){liks[i,7]=1}
				if(x[i]==1 & y[i]==1 & z[i]==1){liks[i,8]=1}
				#If x is ambiguous but the rest are not:
				if(x[i]==2 & y[i]==0 & z[i]==0){liks[i,c(1,2)]=1}
				if(x[i]==2 & y[i]==1 & z[i]==0){liks[i,c(3,5)]=1}
				if(x[i]==2 & y[i]==0 & z[i]==1){liks[i,c(4,6)]=1}
				if(x[i]==2 & y[i]==1 & z[i]==1){liks[i,c(7,8)]=1}
				#If y is ambiguous but the rest are not:
				if(x[i]==0 & y[i]==2 & z[i]==0){liks[i,c(1,3)]=1}
				if(x[i]==1 & y[i]==2 & z[i]==0){liks[i,c(2,5)]=1}
				if(x[i]==0 & y[i]==2 & z[i]==1){liks[i,c(4,7)]=1}
				if(x[i]==1 & y[i]==2 & z[i]==1){liks[i,c(6,8)]=1}
				#If z is ambiguous but the rest are not:
				if(x[i]==0 & y[i]==0 & z[i]==2){liks[i,c(1,4)]=1}
				if(x[i]==0 & y[i]==1 & z[i]==2){liks[i,c(3,7)]=1}
				if(x[i]==1 & y[i]==0 & z[i]==2){liks[i,c(2,6)]=1}
				if(x[i]==1 & y[i]==1 & z[i]==2){liks[i,c(5,8)]=1}
				#If x and y is ambiguous but z is not:
				if(x[i]==2 & y[i]==2 & z[i]==0){liks[i,c(1,2,3,5)]=1}
				if(x[i]==2 & y[i]==2 & z[i]==1){liks[i,c(4,6,7,8)]=1}
				#If x and z is ambiguous but y is not:
				if(x[i]==2 & y[i]==0 & z[i]==2){liks[i,c(1,2,4,6)]=1}
				if(x[i]==2 & y[i]==1 & z[i]==2){liks[i,c(3,5,7,8)]=1}
				#If y and z is ambiguous but x is not:
				if(x[i]==0 & y[i]==2 & z[i]==2){liks[i,c(1,3,4,7)]=1}
				if(x[i]==1 & y[i]==2 & z[i]==2){liks[i,c(2,5,6,8)]=1}
				#All states are ambiguous:			
				if(x[i]==2 & y[i]==2 & z[i]==2){liks[i,1:8]=1}
			}
		}
		Q <- matrix(0, nl^k, nl^k)
		tranQ <- matrix(0, nl^k, nl^k)
	}
	Q[] <- c(p, 0)[rate]
	diag(Q) <- -rowSums(Q)
	phy <- reorder(phy, "pruningwise")
	TIPS <- 1:nb.tip
	anc <- unique(phy$edge[,1])
	comp.branch <- as.list(numeric(nb.tip + nb.node))
	for(i in 1:(nb.tip + nb.node)){
		comp.branch[[i]] <- matrix(0, 1, (dim(liks)[2])+1)
	}
	split.interval = 1
	lik.states <- numeric(nb.tip + nb.node)
	for (i in seq(from = 1, length.out = nb.node)) {
		#The ancestral node at row i is called focal:
		focal <- anc[i]
		#Get descendant information of focal:
		desRows<-which(phy$edge[,1]==focal)
		#Get node information for each descendant:
		desNodes<-phy$edge[desRows,2]
		#Initiates a loop to check if any edges are tip edges: 
		for (desIndex in sequence(length(desRows))){
			#If a tip calculate C_y(i) for the tips and stores in liks matrix:
			if(any(desNodes[desIndex]==phy$edge[,1])==FALSE){
				#This is here to initiate things at the tip, using the branch interval as the unit of time:
				if(hrm == TRUE){
					v <- c(rep(1, k*rate.cat))
				}
				if(hrm == FALSE){
					v <- c(rep(1, nl^k))
				}
				branch.interval <- phy$edge.length[desRows[desIndex]]/split.interval
				Pij <- expm(Q * branch.interval, method=c("Ward77"))
				v <- v * liks[desNodes[desIndex],]
				comp.branch[[desNodes[desIndex]]][1,1] <- branch.interval
				for(i in 1:dim(Pij)[1]){
					L <- Pij[i,] * v
					liks[desNodes[desIndex],i] <- max(L)
					comp.branch[[desNodes[desIndex]]][1,i+1] <- which(L==max(L))[1]
				}	
				#Now that we have the likelihood, we can proceed down the branch collecting the likeliest states:
				if((phy$edge.length[desRows[desIndex]] - branch.interval) > 0){
					branch.calc <- BranchCalc(Q=Q, branch.length=phy$edge.length[desRows[desIndex]]-branch.interval, branch.interval=branch.interval, liks=liks[desNodes[desIndex],], comp.branch=comp.branch[[desNodes[desIndex]]])
					comp.branch[[desNodes[desIndex]]] <- branch.calc$comp.branch
					liks[desNodes[desIndex],] <- branch.calc$liks.final
				}
			}else{
				#Not a tip, but an internal node, and since we already used an interval to calculate things at the node (see below), the branch length is shortened to account for this:
				branch.interval <- phy$edge.length[desRows[desIndex]]/split.interval
				if((phy$edge.length[desRows[desIndex]] - branch.interval) > 0){
					branch.calc <- BranchCalc(Q=Q, branch.length=phy$edge.length[desRows[desIndex]]-branch.interval, branch.interval=branch.interval, liks=liks[desNodes[desIndex],], comp.branch=comp.branch[[desNodes[desIndex]]])				
					comp.branch[[desNodes[desIndex]]] <- branch.calc$comp.branch
					liks[desNodes[desIndex],] <- branch.calc$liks.final
				}
			}
		}
		branch.length <- phy$edge.length[which(phy$edge[,2] == focal)]	
		if(length(branch.length)==0){
			#The focal node is the root, calculate P_k:
			root.state=1
			for (desIndex in sequence(length(desRows))){
				#This is the basic marginal calculation:
				root.state <- root.state * liks[desNodes[desIndex],]
			}
			equil.root <- NULL
			for(i in 1:ncol(Q)){
				posrows <- which(Q[,i] >= 0)
				rowsum <- sum(Q[posrows,i])
				poscols <- which(Q[i,] >= 0)
				colsum <- sum(Q[i,poscols])
				equil.root <- c(equil.root,rowsum/(rowsum+colsum))
			}
			if (is.null(root.p)){
				flat.root = equil.root
				k.rates <- 1/length(which(!is.na(equil.root)))
				flat.root[!is.na(flat.root)] = k.rates
				flat.root[is.na(flat.root)] = 0
				liks[focal, ] <- root.state
			}
			else{
				if(is.character(root.p)){
					# root.p==yang will fix root probabilities based on the inferred rates: q10/(q01+q10), q01/(q01+q10), etc.
					if(root.p == "yang"){
						Q.tmp <- Q
						diag(Q.tmp) <- 0
						root.p <- colSums(Q.tmp) / sum(Q.tmp)
						liks[focal, ] <- root.state * root.p
					}else{
						# root.p==maddfitz will fix root probabilities according to FitzJohn et al 2009 Eq. 10:
						root.p <- root.state / sum(root.state)
						liks[focal,] <- root.state * root.p
					}
				}
				# root.p!==NULL will fix root probabilities based on user supplied vector:
				else{
					liks[focal, ] <- root.state * root.p
				}
			}
		}
		#All other internal nodes, except the root:
		else{
			#Calculates P_ij(t_z):
			branch.interval <- branch.length/split.interval
			Pij <- expm(Q * branch.interval, method=c("Ward77"))
			#Calculates L_z(i):
			if(hrm == TRUE){
				v <- c(rep(1, k*rate.cat))
			}
			if(hrm == FALSE){
				v <- c(rep(1, nl^k))
			}
			for (desIndex in sequence(length(desRows))){
				v <- v * liks[desNodes[desIndex],]
			}
			#Finishes L_z(i):
			comp.branch[[focal]][1,1] <- branch.interval
			for(i in 1:dim(Pij)[1]){
				L <- Pij[i,] * v
				liks[focal,i] <- max(L)
				comp.branch[[focal]][1,i+1] <- which(L==max(L))[1]
			}									
			if(sum(liks[focal,]) < 1e-100){
				#Kicks in arbitrary precision calculations: 
				#library(Rmpfr)
				liks <- mpfr(liks, 15)
			}
		}
	}
	root <- nb.tip + 1L
	comp.branch[[root]][,2:(dim(liks)[2]+1)] <- which(liks[root,]==max(liks[root,]))
	logl <- as.numeric(log(liks[root,which(liks[root,]==max(liks[root,]))]))
	if(include.maps==TRUE){
		obj$loglik <- logl
		obj$mapped.tree<-MapSum(phy, data.sort, liks, comp.branch)
		return(obj)
	}else{
		return(logl)
	}
}


BranchCalc <- function(Q, branch.length, branch.interval, liks, comp.branch){
	comp.tmp <- comp.branch
	liks.tmp <- matrix(liks, 1, 2)
	time.left <- branch.length
	if(time.left == 0){
		v <- c(rep(1, 2))
		Pij <- expm(Q * branch.length, method=c("Ward77"))
		v = v * liks.tmp
		comp.tmp[1,1] <- branch.interval
		for(i in 1:dim(Pij)[1]){
			L <- Pij[i,] * v
			liks.tmp[1,i] <- max(L)
			comp.tmp[1,i+1] <- which(L==max(L))[1]
		}
		comp.branch <- rbind(comp.branch, comp.tmp)
	}else{
		branch.breaks <- ceiling(branch.length / branch.interval)
		for(i in 1:branch.breaks){
			v <- c(rep(1, 2))
			Pij <- expm(Q * branch.interval, method=c("Ward77"))
			v = v * liks.tmp
			comp.tmp[1,1] <- branch.interval
			for(i in 1:dim(Pij)[1]){
				L <- Pij[i,] * v
				liks.tmp[1,i] <- max(L)
				comp.tmp[1,i+1] <- which(L==max(L))[1]
			}
			time.left = time.left - branch.interval
			comp.branch <- rbind(comp.branch, comp.tmp)
			if(time.left < branch.interval){
				branch.interval=time.left
			}
		}
	}
	obj <- NULL
	obj$liks.final <- liks.tmp
	obj$comp.branch <- comp.branch
	return(obj)
}


MapSum <- function(phy, data.sort, liks, comp.branch){
	phy <- reorder(phy, "cladewise")
	mtree<-phy 
	nb.tip = Ntip(mtree)
	nb.node = Nnode(mtree)
	mtree$maps <- as.list(numeric(nrow(phy$edge)))
	mtree$mapped.edge<-matrix(0,nrow(phy$edge),ncol(liks),dimnames=list(paste(phy$edge[,1],",",phy$edge[,2],sep=""),levels(as.factor(data.sort[,1]))))
	#A throw-away vector to keep track of end states of where we have been:
	comp <- numeric(nb.tip + nb.node)
	root <- nb.tip + 1L
	comp[root] <- comp.branch[[root]][1,2]
	#reverse the order of the focal nodes:
	anc <- unique(phy$edge[,1])
	count <- 0
	for (i in seq(from = 1, length.out = nb.node)){
		#The ancestral node at row i is called focal:
		focal <- anc[i]
		parent.state <- comp[focal]
		#Get descendant information of focal:
		desRows <- which(phy$edge[,1]==focal)
		#Get node information for each descendant:
		desNodes<-phy$edge[desRows,2]
		for(desIndex in sequence(length(desRows))){
			count = count+1
			previous.state = parent.state
			map <- 0
			names(map) <- parent.state
			sampled.points <- dim(comp.branch[[desNodes[desIndex]]])[1]
			map.count <- 1
			for(k in sampled.points:1){
				current.state <- comp.branch[[desNodes[desIndex]]][k,previous.state+1]
				if(current.state == previous.state){
					map[map.count] <- map[map.count] + comp.branch[[desNodes[desIndex]]][k,1]
				}else{
					new.map <- comp.branch[[desNodes[desIndex]]][k,1]
					names(new.map) <- current.state
					map <- c(map, new.map)
					previous.state <- current.state
					map.count <- map.count + 1
				}
				if(k == 1){
					comp[desNodes[desIndex]] <- previous.state
				}				
			}
			proper.spot <- which(desNodes[desIndex]==phy$edge[,2])
			#Maps need to be in the order of the edge matrix:
			mtree$maps[[proper.spot]]<-map
		}
	}
	#Taken from the make.simmap function in a recent version of phytools to ensure compatibility with plotSimmap:
	allstates<-vector()
	for(j in 1:nrow(mtree$edge)) allstates<-c(allstates,names(mtree$maps[[j]]))
	allstates<-unique(allstates)
	mtree$mapped.edge<-matrix(data=0,length(mtree$edge.length),length(allstates),dimnames=list(apply(mtree$edge,1,function(x) paste(x,collapse=",")),state=allstates))
	for(j in 1:length(mtree$maps)) for(k in 1:length(mtree$maps[[j]])) mtree$mapped.edge[j,names(mtree$maps[[j]])[k]]<-mtree$mapped.edge[j,names(mtree$maps[[j]])[k]]+mtree$maps[[j]][k]
	return(mtree)
}


#########################
#    match.tree.data    #
#########################
# Compares a tree and data to make sure they include the same taxa
# Taxa which are in the tree, but not the data matrix, are added to the matrix and coded as missing data.
# Any taxa in the data matrix which are not in the tree are removed from the matrix
# The function returns an object with three parts:
#	$phy: the tree
#	$data: the matrix, omitting taxa not in tree and taxa that were present in the tree but not in the matrix
#	$message.data: a brief message explaining modifications (if any) to the data
#	$message.tree: a brief message explaining modificatoins (if any) to the tree
match.tree.data <- function(phy, data){
	matchobj <- NULL
	matchobj$phy <- phy
	matchobj$data <- data
	matchobj$message.data <- NULL
	matchobj$message.tree <- NULL
	# First look at data matrix to see if each taxon in matrix is also in tree
	missing.fromtree <- NULL
	for(datarow in 1:length(data[,1])){
		if(is.na(match(data[datarow,1],phy$tip.label))){
			missing.fromtree <- c(missing.fromtree,datarow)
		}
	}
	if(length(missing.fromtree) > 0){ # At least one taxa is listed in the matrix, but is not in the tree
		# Make message so user knows taxa have been removed
		matchobj$message.data <- "The following taxa in the data matrix were not in the tree and were excluded from analysis: "
		first <- TRUE
		for(toRemove in 1:length(missing.fromtree)){
			if(first){
				matchobj$message.data <- paste(matchobj$message.data,as.character(data[missing.fromtree[toRemove],1]),sep="")
				first <- FALSE
			} else { #not the first one, so add leading comma
				matchobj$message.data <- paste(matchobj$message.data,", ",as.character(data[missing.fromtree[toRemove],1]),sep="")
			}
		}
		matchobj$data <- data[-missing.fromtree,] # omits those data rows which have no match in the tree
		for(datacol in 2:length(matchobj$data[1,])){
			matchobj$data[,datacol] <- factor(matchobj$data[,datacol]) # have to use factor to remove any factors not present in the final dataset
		}
	}
	
	missing.taxa <- NULL
	for(tip in 1:length(phy$tip.label)){
		if(is.na(match(phy$tip.label[tip],matchobj$data[,1]))){
			if(is.null(matchobj$message.tree)){ # The first missing taxon
				missing.taxa <- as.character(phy$tip.label[tip])
				matchobj$message.tree <- "The following taxa were in the tree but did not have corresponding data in the data matrix.  They are coded as missing data for subsequent analyses: "
			} else { # not the first missing taxon, add with leading comma
				missing.taxa <- paste(missing.taxa,", ",as.character(phy$tip.label[tip]),sep="")
			}
			# missing taxa will be coded as having missing data "?"
			addtaxon <- as.character(phy$tip.label[tip])
			numcols <- length(matchobj$data[1,])
			newrow <- matrix(as.character("\x3F"),1,numcols) # absurd, but it works
			newrow[1,1] <- addtaxon
			newrowdf <- data.frame(newrow)
			colnames(newrowdf) <- colnames(matchobj$data)
			matchobj$data <- rbind(matchobj$data,newrowdf)
		}
	}
	rownames(matchobj$data) <- matchobj$data[,1] # Use first column (taxon names) as row names
	matchobj$data <- matchobj$data[matchobj$phy$tip.label,] # Sort by order in tree
	rownames(matchobj$data) <- NULL # remove row names after sorting
	if(!is.null(missing.taxa)){
		matchobj$message.tree <- paste(matchobj$message.tree,missing.taxa,sep="")
	}
	return(matchobj)
}


##############
#  findAmps  #
##############
# A function to find positions of ampersands for separating different states.  
# Will allow character state to be greater than one character long.
findAmps <- function(string, charnum){
	if(!is.character(string)) return(NULL)
	locs <- NULL # Will hold location values
	for(charnum in 1:nchar(as.character(string))){
		if(substr(string,charnum,charnum) == "&"){
			locs <- c(locs,charnum)
		}
	}
	return(locs)
}

##############
# factorData #
##############
# Function to make factored matrix as levels are discovered.
factorData <- function(data,whichchar=1,charnum){
	charcol <- whichchar+1
	factored <- NULL # will become the matrix.  Starts with no data.
	lvls <- NULL
	numrows <- length(data[,charcol])
	missing <- NULL
	
	for(row in 1:numrows){
		currlvl <- NULL
		levelstring <- as.character(data[row,charcol])
		ampLocs <- findAmps(levelstring, charnum)
		if(length(ampLocs) == 0){ #No ampersands, character is monomorphic
			currlvl <- levelstring
			if(currlvl == "?" || currlvl == "-"){ # Check for missing data
				missing <- c(missing,row) # add to list of taxa with missing values, will fill in entire row later
			}
			else { # Not missing data
				if(length(which(lvls == currlvl)) == 0){# encountered a level not seen yet
					if(length(factored) == 0){ # Matrix is empty, need to create it
						factored <- matrix(0,numrows,1)
						colnames(factored) <- currlvl
						rownames(factored) <- rownames(data)
					} else { # matrix already exists, but need to add a column for the new level
						zerocolumn <- rep(0,numrows)
						factored <- cbind(factored, zerocolumn)
						colnames(factored)[length(factored[1,])] <- currlvl
					}
					lvls <- c(lvls,currlvl) # add that level to the list
				} # already found this level in another state.  Set the value to one
				whichlvl <- which(lvls == currlvl) # this index number should correspond to the column number of the state
				factored[row,whichlvl] <- 1
			}
		} else { #At least one ampersand found, polymorphic character
			start <- 1
			numlvls <- length(ampLocs)+1
			for(part in 1:numlvls){
				# Pull out level from levelstring
				if(part <= length(ampLocs)){ # Havent reached the last state
					currlvl <- substr(levelstring,start,(ampLocs[part]-1)) # pull out value between start and the location-1 of the next ampersand
				} else { # Final state in list
					currlvl <- substr(levelstring,start,nchar(levelstring)) # pull out value between start and the last character of the string
				}
				if(currlvl == "?" || currlvl == "-"){ # Missing data, but polymorphic?
					missing <- c(missing,row) # add to list of taxa with missing values, will fill in entire row later
				}
				else { # Not missing data
					if(length(which(lvls == currlvl)) == 0){# encountered a level not seen yet
						if(length(factored) == 0){ # Matrix is empty, need to create it
							factored <- matrix(0,numrows,1)
							colnames(factored) <- currlvl
							rownames(factored) <- rownames(data)
						} else { # matrix already exists, but need to add a column for the new level
							zerocolumn <- rep(0,numrows)
							factored <- cbind(factored, zerocolumn)
							colnames(factored)[length(factored[1,])] <- currlvl
						}
						lvls <- c(lvls,currlvl) # add that level to the list
					} # already found this level in another state.  Set the value to one
					whichlvl <- which(lvls == currlvl) # this index number should correspond to the column number of the state
					factored[row,whichlvl] <- 1
					start <- ampLocs[part] + 1
				}
			}
		}
	}
	#Need to deal with any rows with missing data; fill in NA for all columns for that row
	for(missingrows in 1:length(missing)){
		for(column in 1:length(factored[1,])){
			factored[missing[missingrows],column] <- 1 # All states equally likely
		}
	}
	factored <- factored[,order(colnames(factored))]
	return(factored)
}


#library(phytools)
#library(corHMM)
#library(expm)
#Pr.mat <- matrix(c(0.70, 0.45, 0.30, 0.55), 2, 2)
#Q.mat <- logm(Pr.mat)
#rates <- c(Q.mat[2,1], Q.mat[1,2])
#phy<-read.tree("pupko.tre")
#trait<-read.delim("pupko.data.txt")
#pp <- PathLik(phy, trait, p=rates, ntraits=1, charnum=1, model="ARD", root.p="yang")



