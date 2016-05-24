"pd.calc" <-
function(cm,  tip.subset=NULL, method="TBL", root.edge=FALSE){

	# check we have a valid method
	method <- match.arg(method, c("TBL", "MST", "UEH","SBL", "TIP"))
	
	# check we have a clade matrix and, if not, get one
	if(class(cm) != "clade.matrix"){
		if(class(cm) == "phylo"){
			warning("Converting phylo object to clade.matrix object")
			cm <- clade.matrix(cm)
		} else { 
			stop("pd.calc requires a phylogeny")	
		}
	}

    # if requested, drop the root edge
    nSpp <- dim(cm$clade.matrix)[2]
    if(! root.edge) cm$edge.length[nSpp + 1] <- 0
    	
	# subset the tips if requested
	if(!is.null(tip.subset)){
		# could be names - in which case they must match tip labels
		# could be numbers - in which case they must be in range
		switch(mode(tip.subset),
			"character" = {
				tip.subset <- match(tip.subset, cm$tip.label)
				if(any(is.na(tip.subset))){
					stop("Unmatched names in tip.subset")}},
			"numeric" = {
				if(any(tip.subset %in% 1:dim(cm$clade.matrix)[2] == FALSE)){
					stop("numeric tip.subset contains outside the range 1 to number of tips")}},
			stop("tip.subset must be either a vector of names or numbers"))
	} else {
		tip.subset <- 1:dim(cm$clade.matrix)[2]
	}
	
	#choose method
	switch(method,
		"TBL" = {
			edge.in.matrix <- cm$clade.matrix[,tip.subset]
			if(is.array(edge.in.matrix)){
				edge.in.matrix <- rowSums(edge.in.matrix)}
			edge.in.matrix <- edge.in.matrix > 0 
		},
		"MST" = {
			edge.in.matrix <- cm$clade.matrix[, tip.subset]
			if(is.array(edge.in.matrix)){
				edge.in.matrix <- rowSums(edge.in.matrix)}
			edge.in.matrix <- edge.in.matrix > 0 & edge.in.matrix < length(tip.subset)
		},
		"TIP" = {
			edge.in.matrix <- cm$clade.matrix[,tip.subset]
			if(is.array(edge.in.matrix)){
				edge.in.matrix <- rowSums(edge.in.matrix)}
			edge.in.matrix <- edge.in.matrix > 0 
			edge.in.matrix[(dim(cm$clade.matrix)[2]+1):length(edge.in.matrix)] <- FALSE
		},
		"UEH" = {
			edge.in.matrix <- cm$clade.matrix[,tip.subset]
			if(is.array(edge.in.matrix)){
				edge.in.matrix <- rowSums(edge.in.matrix)}
			edge.in.matrix <- edge.in.matrix == 1 
		},
        "SBL" = {
			edge.in.matrix <- cm$clade.matrix[,tip.subset]
			if(is.array(edge.in.matrix)){
				edge.in.matrix <- rowSums(edge.in.matrix)}
			edge.in.matrix <- edge.in.matrix > 1 
		})
	
	pd <- sum(cm$edge.len[edge.in.matrix])
	RET <- structure(.Data=pd, pd.method=method)
		
	return(RET)
	
}


"pd.bootstrap" <-
function(cm, ntips, reps=1000, method="TBL", tip.weights=NULL){
    
	# check we have a valid method
	method <- match.arg(method, c("TBL", "MST", "UEH", "SBL", "TIP"))

	# check we have a clade matrix and, if not, get one
	if(class(cm) != "clade.matrix"){
		if(class(cm) == "phylo"){
			warning("Converting phylo object to clade.matrix object")
			cm <- clade.matrix(cm)
		} else { 
			stop("pd.calc requires a phylogeny")	
		}
	}
	
	# check for sensible sample
	total.nb.tips <- dim(cm$clade.matrix)[2]
	if(!(ntips %in% 1:(total.nb.tips - 1))){
		stop("'sample' must be a positive integer lower than the number of tips")}
	
	#set up the store
	pd.store <- numeric(reps)
	tips <- 1:total.nb.tips
	
	# if there are weights make sure they go in the right place...
	if( ! is.null(tip.weights)){
	        
        # if the vector is named then match the order to tips
        if(! is.null(names(tip.weights))) {
            wght.match <- match(cm$tip.label, names(tip.weights))
            
            if(any(is.na(wght.match))){ # this is not elegant but can't work out how to stop and return
                warning("The returned tip labels have no matching named element in tip.weights")
                return(cm$tip.label[is.na(wght.match)])
                }
            
            tip.weights <- tip.weights[wght.match]
            
        } else {
             stop("'weights' must be a vector of weights, named to match the tip labels")
        }
    }
    
	# get the pd values
    for(rep in seq(along=pd.store)){
        which.tips <- sample(tips, ntips, prob=tip.weights)
        pd.store[rep] <- pd.calc(cm, tip.subset=which.tips, method=method)
	}
	
	return(structure(.Data=pd.store, pd.method=method))
}

ed.calc <- function(cm, polytomy.cf=c("isaac","mooers","none")){

	#Nick Isaac, March 2009 + David Orme 2011
	#takes the phylogeny and returns a list containing ED scores of a) species and b) branches
	#the polytomy.cf argument specifies which set of polytomies should be applied. There are three options: 
		#"isaac" is as the EDGE paper of 2007, based on logarithmic decay with node size (too harsh on large nodes)
		#"mooers" is empirical, based on a pure-birth process
		#"none" : no correction
	
	#CHANGES TO MAKE:
		#1) Optional data frame containing IUCN categories
		#2) Optional data frame containing richness weights - to account for missing species
	
	# check we have a clade matrix and, if not, get one
	
	if(inherits(cm,"phylo")){
		warning("Converting phylo object to clade.matrix object")
		cm <- clade.matrix(cm)
	} else if(! inherits(cm, 'clade.matrix')) { 
		stop("pd.calc requires a phylogeny")	
	}
	
	polytomy.cf <- match.arg(polytomy.cf)
	
	## get raw edge scores 
	## (allocates equal proportions of each branch length between all descendent species)
	branch <- data.frame(len=cm$edge.length, nSp=rowSums(cm$clade.matrix))
	branch$ED <- with(branch, len/nSp)
	
	## polytomy corrections (the branch lengths at soft polytomies
	## overestimate the amount of evolution going on)
	
	## get the group size at the parent nodes (i.e. the number of siblings at a node)
	branch$parent <- cm$edge[,1][match(rownames(branch), cm$edge[,2])]
	node.size <- as.data.frame(table(cm$edge[,1]))
	branch$node.size <- node.size$Freq[match(branch$parent, node.size$Var1)]
	
	branch$ED.cor <- switch(polytomy.cf,
					'isaac'  = {with(branch,  ifelse(node.size > 57, 0, ED * (1.081 - 0.267 * log(node.size))))},
					'mooers' = {with(branch,  ED/node.size * (node.size - 1)/sapply(node.size,function(n){if(is.na(n)) NA else sum(1/2:n)}))},
					'none'   = branch$ED)
	branch$ED.cor <- with(branch, ifelse(node.size > 2, ED.cor, ED))
	
 	# get species edge sums
	edge.matrix <- branch$ED.cor * cm$clade.matrix
	spp.ED <- data.frame(species=cm$tip.label, ED=colSums(edge.matrix, na.rm=TRUE), stringsAsFactors=FALSE)
	
	return(list(spp=spp.ED,branch=branch))
}


