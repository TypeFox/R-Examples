
comparative.data <- function(phy, data, names.col, vcv=FALSE, vcv.dim=2, na.omit=TRUE, force.root=FALSE, warn.dropped=FALSE, scope=NULL){

    # TODO - is something odd happening with missing arguments?
    
    # record call and dataset names
    phy.name <- deparse(substitute(phy))
    data.name <- deparse(substitute(data))
    names.col <- as.character(substitute(names.col))

    # check inputs are what they should be
    # DATA:
    
        # ...is a dataframe
        if(! is.data.frame(data)) stop("'data' must be an object of class 'data.frame'.")
        # ...contains the name column and make sure it is of mode character
        namesInd <- match(names.col, names(data))
        if(is.na(namesInd)) {
            stop("Names column '",  names.col, "' not found in data frame '", data.name, "'")
        }
        rownames(data) <- as.character(data[,namesInd])
        # drop the names column
        data <- data[,-namesInd, drop=FALSE]
        
    # PHYLOGENY:
        # check the phylogeny is a rooted phylogeny and test for stupid tip labels
        if(! inherits(phy, "phylo")) 
            stop("'", deparse(substitute(phy)), "' not of class 'phylo'")
        if(! is.rooted(phy)){
			if(force.root){
				phy$root.edge <- 1
			} else {
				stop("'", deparse(substitute(phy)), "' is not rooted or has a basal polytomy.")
			}
		}
            
        if(any(duplicated(phy$tip.label))) stop('Duplicate tip labels present in phylogeny')
        if(any(duplicated(c(phy$tip.label, phy$node.label)))) stop('Labels duplicated between tips and nodes in phylogeny')


    # MERGE
    
        # store original dataset size
        origTips <- with(phy, max(edge) - Nnode)
        origData <- nrow(data)
        
        # find the intersection between tip.labels and names in data frame
        in.both <- intersect(rownames(data), phy$tip.label)
        if(length(in.both) == 0 ) stop("No tips are common to the dataset and phylogeny")
        
        # i >> ditch rows with no tip
        row.in.tree <- match(rownames(data), in.both)
        row.not.in.tree <- rownames(data)[is.na(row.in.tree)]
        data <- subset(data, !is.na(row.in.tree))
    
        # ii >> ditch tips which have no rows.
        tip.in.data <-  match(phy$tip.label, in.both)
        to.drop <- phy$tip.label[is.na(tip.in.data)]
        
        #  get subset of phylogeny to be used
        if(length(to.drop) > 0) matchedPhy <- drop.tip(phy, to.drop) else matchedPhy <- phy
        
        # useful info...
        root <- length(matchedPhy$tip.label) + 1

        # get the data into the same order as the tips
        tip.order <- match(matchedPhy$tip.label, rownames(data))
        if(any(is.na(tip.order))) stop("Problem with sorting data frame: mismatch between tip labels and data frame labels")
        data <- data[tip.order,, drop=FALSE]
        
        # Label the data frame rows by tip label
        rownames(data) <- matchedPhy$tip.label
        
		# add or supplement node labels - used as references in contrast calculation
		IntNd <- root:max(matchedPhy$edge)
		if(is.null(matchedPhy$node.label)){
			matchedPhy$node.label <- IntNd
		} else {
			# set up missing node labels
			matchedPhy$node.label <- ifelse(matchedPhy$node.label == "", NA, matchedPhy$node.label)
			if(any(duplicated(na.omit(matchedPhy$node.label)))) stop('Duplicate node labels present in phylogeny')
			matchedPhy$node.label <- ifelse(is.na(matchedPhy$node.label),  IntNd, matchedPhy$node.label)
		}
		
		
    # Compile comparative dataset
	# COULD HAVE phy components as first level slots rather than nested within $phy
	# and have comparative.data inherit methods from phylo, but that presupposes
	# that none of the phylo functions you might use strip out contents
	
    RET <- list(phy=matchedPhy, data = data, 
                data.name=data.name, phy.name=phy.name, 
                dropped=list(tips=to.drop, unmatched.rows=row.not.in.tree))
    class(RET) <- 'comparative.data'
    
    # Add a VCV array if requested
    if(vcv) {
        RET$vcv <- VCV.array(matchedPhy, dim=vcv.dim)
        RET$vcv.dim <- vcv.dim
    }
    
    # NA handling
    if(na.omit){
    	before.drop.rows <- rownames(RET$data)
        RET <- na.omit(RET, scope)
        if(!identical(rownames(RET$data), before.drop.rows)) RET$dropped$NA.rows <- before.drop.rows
    }
    
	if(warn.dropped){
		if(any(sapply(RET$dropped, length) > 0)) warning('Data dropped in compiling comparative data object')
	}

    return(RET)
}

# some useful generics
print.comparative.data <- function(x, ...){

    # basic summary data
    cat("Comparative dataset of", nrow(x$data), "taxa:\n")
	if(! is.null(attr(x, 'growTree')))     cat("Simulated using growTree: object contains node data and simulation details\n")
    cat("Phylogeny:", x$phy.name, "\n")
    cat("   ", length(x$phy$tip.label), " tips, ", x$phy$Nnode, " internal nodes\n  ", sep='')
    # this is a bit of a hack - can't get str for a vector to take an indent
    str(x$phy$tip.label)
    if(! is.null(x$vcv)){
	    cat('VCV matrix present:\n  ')
	    str(x$vcv, give.attr=FALSE)
	}
    cat("Data:" , x$data.name, "\n")
    str(as.list(x$data), no.list=TRUE, indent='   ')

	# report on mismatch on merge
	dropCount <- sapply(x$dropped, length)
    if(any(dropCount)){
	    cat('Dropped taxa:\n')
	    cat('   ', x$phy.name , ' { ', dropCount[1], ' ( ',nrow(x$data), 
	        ' } ', dropCount[2], ' ) ', x$data.name, sep='')
    }
	if(! is.null(attr(x, 'na.omit.scope'))){
		cat('\nScope of complete data:\n', deparse(attr(x, 'na.omit.scope')))
	}
}

na.omit.comparative.data <- function(object, scope=NULL, ...){

    # strips data rows, tips and vcv row/columns for a comparative.data object
	if(! is.null(scope)){
		if(! is.null(attr(object, 'na.omit.scope'))) stop('Scope of comparative data set already set.')
		if(! inherits(scope, 'formula')) stop('scope must be a model formula.')
		vars <- all.vars(scope)
		if(any(is.na(match(vars, names(object$data))))) stop('Variables in scope not in data.')
		to.drop <- which(! complete.cases(object$data[, vars]))
		attr(object, 'na.omit.scope') <- scope
	} else {
    	to.drop <- which(! complete.cases(object$data))
	}
	# test for completely empty dataset
	if(length(to.drop) == nrow(object$data)) warning('No complete rows present in data.')
	
    # guard against ape 'feature' of dropping all tips for empty vectors
    if(length(to.drop) > 0){
        
        # lose bits of tree
        object$phy <- drop.tip(object$phy, to.drop) 
        # lose rows
        object$data <- object$data[-to.drop,, drop=FALSE]
        
        # lose VCV elements if needed
        if(! is.null(object$vcv)){ 
			if((object$vcv.dim == 2)){
				object$vcv <- object$vcv[-to.drop, -to.drop] 
			} else {
				object$vcv <- object$vcv[-to.drop, -to.drop, ] 
			}
		}
    
		# add to dropped list
		object$dropped$tips <- c(object$dropped$tips, to.drop)
	}
    
	if(! is.null(attr(object, 'growTree'))) cat("Warning: subsetting of node data not implemented.\n")
    return(object)
}

subset.comparative.data <- function(x, subset, select,  ...){

    ## ripping out the innards of subset.data.frame
    if (missing(subset)) 
        r <- TRUE
    else {
        e <- substitute(subset)
        r <- eval(e, x$data, parent.frame())
        if (!is.logical(r)) 
            stop("'subset' must evaluate to logical")
        r <- r & !is.na(r)
    }
    if (missing(select)) 
        vars <- TRUE
    else {
        nl <- as.list(seq_along(x$data))
        names(nl) <- names(x$data)
        vars <- eval(substitute(select), nl, parent.frame())
    }
    
    ## now know which rows and columns to keep in the data frame 
    ## guard against ape 'feature' of dropping all tips for empty vectors
    if(any(! r)){
        
        to.drop <- which(! r)
        # lose bits of tree
        x$phy <- drop.tip(x$phy, to.drop) 
        # lose rows
        x$data <- x$data[r, vars, drop = FALSE] # CANNOT lose data.frame-ness

        # lose VCV elements if needed
        if(! is.null(x$vcv)){ 
			if((x$vcv.dim == 2)){
				x$vcv <- x$vcv[r, r] 
			} else {
				x$vcv <- x$vcv[r, r, ] 
			}
		}
		
		# add to dropped list
		x$dropped$tips <- c(x$dropped$tips, to.drop)
    }

	if(! is.null(attr(x, 'growTree'))) cat("Warning: subsetting of node data not implemented.\n")

    return(x)
}

"[.comparative.data" <- function(x, i, j) {
	
	# how many args?
	# 2 and missing(i) = x[] --> return x untouched
	# 2 and ! missing(i) = x[i] --> column subset 
	# otherwise 3 = x[i,j] or x[i,] or x[,j] or x[,] 
	if( nargs() ==  2) {
		if(missing(i)){
			return(x)
		} else {
			j <- i
			hasI <- FALSE
		}
	} else { hasI <- TRUE }
	
	# no drop argument permitted. can't lose dataframeness
	if(! missing(j)){
		if(is.null(j)) stop('Null indices not permitted on comparative data objects')
		x$data <- x$data[,j, drop=FALSE]
	}

	# no recycling, no out of index rows
	# no simple reordering possible because the tree implies
	# an order to the data frame
	if(! missing(i) & hasI){
		if(is.null(i)) stop('Null indices not permitted on comparative data objects')
		rownames <- x$phy$tip.label
		if(is.character(i)){
			toKeep <- na.omit(match(i, rownames))
			if(length(toKeep) != length(i)) warning('Some tip names were not found')
			toKeep <- rownames[toKeep]
		} else if (is.numeric(i)) {
			# convert to integer (same as [.data.frame)
			i <- as.integer(i)
			if(all(i > 0)){
				toKeep <- intersect(i, seq_along(rownames))
			} else if(all(i < 0)) {
				toKeep <- setdiff(seq_along(rownames), abs(i))
			} else {
				stop("only 0's may be mixed with negative subscripts")
			}
			toKeep <- rownames[toKeep]
			if(! all(abs(i) %in% seq_along(rownames))) warning('Some row numbers were not found')
		} else if (is.logical(i)) {
			if(length(i) != length(rownames)) stop('Logical index does not match number of tips')
			toKeep <- rownames[i]
		}
		
		# Work out which to drop
		#  - 'drop.tip' will create a pointless tree if asked to drop nothing from a tree 
		rowToKeep <- match(toKeep, rownames)
		if(length(rowToKeep) < 2) stop('Comparative dataset cannot be reduced to fewer than 2 taxa')
		if(length(rowToKeep) != length(rownames)) {
			toDrop    <- setdiff(rownames, toKeep)
			# this assumes that drop.tip preserves the remaining order - testing suggests ok
			if(length(rowToKeep) == length(rownames)) warning('No taxa to drop from comparative dataset') else x$phy <- drop.tip(x$phy, toDrop)
			
			# add to dropped list
			x$dropped$tips <- c(x$dropped$tips, toDrop)
			
        	# lose VCV elements if needed
        	if(! is.null(x$vcv)){ 
				if((x$vcv.dim == 2)){
					x$vcv <- x$vcv[rowToKeep, rowToKeep] 
					} else {
						x$vcv <- x$vcv[rowToKeep, rowToKeep, ] 
					}
				}
			}
        x$data <- x$data[rowToKeep,, drop=FALSE]
	}
	
	if(! is.null(attr(x, 'growTree'))) cat("Warning: subsetting of node data not implemented.\n")
	return(x)
}

reorder.comparative.data <- function(x, order = "cladewise", ...){
	
	# Uses ape reorder code
	order <- match.arg(order, c("cladewise", "pruningwise"))

	# test for existing order to avoid duplicate calls
    if (!is.null(attr(x, "order"))) 
        if (attr(x, "order") == order) 
            return(x)
	
	# exclude 2 taxon trees
    nb.node <- x$phy$Nnode
    if (nb.node == 1) 
        return(x)
	
	# otherwise
    nb.tip <- length(x$phy$tip.label)
    nb.edge <- dim(x$phy$edge)[1]

	# get the new order - testing for ape version for change in code
	if( packageVersion('ape') >= '3.0.5'){
		neworder <- reorder(x$phy, order=order, index.only=TRUE)
	} else {
		stop('ape versions < 3.0.4 no longer supported.')
	   # neworder <- if (order == "cladewise") 
	   #     .C("neworder_cladewise", as.integer(nb.tip), as.integer(x$phy$edge[, 
	   #         1]), as.integer(x$phy$edge[, 2]), as.integer(nb.edge), 
	   #         integer(nb.edge), PACKAGE = "ape")[[5]]
	   # else .C("neworder_pruningwise", as.integer(nb.tip), as.integer(nb.node), 
	   #     as.integer(x$phy$edge[, 1]), as.integer(x$phy$edge[, 2]), as.integer(nb.edge), 
	   #     integer(nb.edge), PACKAGE = "ape")[[6]]
	}
	
	# apply new order to elements
    x$phy$edge <- x$phy$edge[neworder, ]
    attr(x$phy, "order") <- order

    if (!is.null(x$phy$edge.length)) 
        x$phy$edge.length <- x$phy$edge.length[neworder]

    # lose VCV elements if needed
    if(! is.null(x$vcv)){ 
		if((x$vcv.dim == 2)){
			x$vcv <- x$vcv[neworder, neworder] 
		} else {
			x$vcv <- x$vcv[neworder, neworder, ] 
		}
	}
	
	x$dat <- x$dat[neworder,]
    attr(x, "order") <- order
    
	return(x)
}

caicStyleArgs <- function(phy, data, names.col, warn.dropped=FALSE){
	
	# general function to handle old style non-'comparative.data' use of 
	# crunch, brunch, macrocaic, phylo.d functions 
	
	if(missing(data)) stop('data object is missing')
	if(missing(phy)) stop('phy object is missing')

	# check the classes (and allow for them being in the wrong order)
	args <- list(data, phy)
	argClass <- sapply(args, class)
	
	# bail back to calling function if we don't have targets
	# i.e. precisely a phylogeny and a data.frame
	targets <- c('data.frame', 'phylo')
	if(! identical( sort(intersect( argClass, targets)), targets)){
		return(NULL)
	}
	
	# try and build a comparative data object
	#  - allow for the order of the variables to be reversed: phylo.d, I'm looking at you
	
    if(argClass[1] == 'data.frame'){
		data <- eval(substitute(comparative.data(phy = phy, data = data, names.col = XXX, warn.dropped=warn.dropped), list(XXX=names.col)))
	} else {
		data <- eval(substitute(comparative.data(phy = data, data = phy, names.col = XXX, warn.dropped=warn.dropped), list(XXX=names.col)))
	}
	return(data)
	
}

## ## THIS NEEDS SOME WORK TO PASS ALL THESE
## crunchMod <- crunch(Egg.Mass ~ F.Mass + M.Mass, phy=shorebird.tree, data=shorebird.data, names.col=Species)
## crunchMod <- crunch(Egg.Mass ~ F.Mass + M.Mass, shorebird.tree, shorebird.data, names.col=Species)
## crunchMod <- crunch(Egg.Mass ~ F.Mass + M.Mass, shorebird.data, shorebird.tree, names.col=Species)
## ## BREAKS WITH THE FOLLOWING MISSING ARGUMENTS BUT THESE WOULD HAVE BROKEN ANYWAY
## crunchMod <- crunch(Egg.Mass ~ F.Mass + M.Mass, phy=shorebird.tree, names.col=Species)
## crunchMod <- crunch(Egg.Mass ~ F.Mass + M.Mass, data=shorebird.data, names.col=Species) 
## crunchMod <- crunch(Egg.Mass ~ F.Mass + M.Mass, shorebird.tree, names.col=Species)
## crunchMod <- crunch(Egg.Mass ~ F.Mass + M.Mass, shorebird.data, names.col=Species)


as.comparative.data <- function(x, ...){
    
	if(inherits(x, 'comparative.data')){
		return(x)
	} else {
		UseMethod('as.comparative.data')
	}
}

as.comparative.data.growTree <- function(x, ...){
	
    lineages <- x$lineages
    clade <- x$clade
    
	# get a phylo object
	
    if(dim(lineages)[1] > 1){
        # extract information (excluding root node)
        linNoR <- lineages[-1,]
            
        # now need to coerce the numbering into ape style
        parentMap <- data.frame(linPar=unique(linNoR$parent.id), apePar=with(clade, (nTip+1):nLin))
        childMap <- data.frame(linPar=with(linNoR, id[tip]), apePar=with(clade, 1:nTip))
        nodeMap <- rbind(parentMap, childMap)
        linNoR$pnode <- nodeMap$apePar[match(linNoR$parent.id, nodeMap$linPar)]
        linNoR$node <- nodeMap$apePar[match(linNoR$id, nodeMap$linPar)]
        
        edge <- as.matrix(linNoR[, c("pnode", "node")])
        dimnames(edge) <- NULL
        
        # lets be honest... the caic.code column is really only there as a cheap
        # route to a cladewise sorting of the edge matrix!
        ord <- order(linNoR$caic.code)
        edge <- edge[ord,]
        edge.length <- linNoR$lin.age[ord]
        
        phy <- list(edge=edge, edge.length=edge.length, tip.label=1:clade$nTip, 
                    Nnode=with(clade, nLin-nTip), root.edge=lineages$lin.age[1])              
        class(phy) <- "phylo"
        
    } else {
        # have an unspeciated root so put in slightly  as a single tip with a root edge of zero
        phy <- list(edge=matrix(c(1,2), ncol=2), edge.length=lineages$lin.age[1], 
                    tip.label=1, root.edge=0, Nnode=1)
        class(phy) <- "phylo"
    }
    
	# sort out comparative data
	lastRules <- x$epochRules[[length(x$epochRules)]]
	datCol <- c('node', 'lin.age', 'birth.time', 'death.time', 'extinct', 'tip',
				names(lastRules$ct.set$ct.start), names(lastRules$dt.rates))
	dat <- linNoR[,datCol]
    
	# get tip data set without tips flag
	tipData <- dat[dat$tip, -6]
	nodeData <- dat[! dat$tip, -5:-6]
 	com <- comparative.data(phy, tipData, 'node', na.omit=FALSE)
	com$node.data <- nodeData
    com$epochRules <- x$epochRules
    attr(com, 'growTree') <- TRUE

    return(com)
}

## as.comparative.data.phylo4d <- function(x, ...){
## 	
## 	
## 	
## }

## x <- comparative.data(shorebird.tree, shorebird.data, 'Species')
## x[]
## x[,]
## x[2:3]
## x[, 2:3]
## x[1:15, ]
## x[1:15, 2:3]

## ## $ method: Don't think this is possible with S3 - want $ to be able
## ## to give back a column from data - without breaking the use
## ## of $ for subsetting. Could use name matching to figure out which
## ## but then duplicate names in data and in the class list are a huge
## ## programming gotcha.

## "$.comparative.data" <- function(x, name) {
## 	
## 	# careful to avoid using $ in here otherwise infinite recursion kicks off.
## 	return(x[['data']][, name])
## 
## }
