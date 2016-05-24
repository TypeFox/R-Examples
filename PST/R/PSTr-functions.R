
## If node not in keep.list it is tagged as pruned
node.keep <- function(x, keep.list, clist) {
	## 
	if (!is.null(clist)) {
		tmp <- unlist(lapply(clist, node.parent))
	}	

	if (!x@path %in% keep.list & (x@leaf | (!is.null(clist) && !x@path %in% tmp))) {
		x@pruned[] <- TRUE
	}
	return(x)
}

## Return the label of a parent's node 
node.parent <- function(x, segmented=FALSE) {
	if (x@order==1) { 
		parent <- "e" 
	} else {
		parent <- unlist(strsplit(x@path,"-"))
		parent <- parent[2:length(parent)]
		parent <- paste(parent, collapse="-")
	}
	
	if (segmented) {
		parent <- cbind(parent, group=rownames(x@prob))
	}

	return(parent)
}

## gain function
node.gain <- function(x, plist, gain, C, clist) {
	parent <- plist[[node.parent(x)]]

	## 
	if (!is.null(clist)) {
		tmp <- unlist(lapply(clist, node.parent))
		children <- clist[which(tmp==x@path)]
		if (length(children)>0) {
			children <- unlist(lapply(children, function(x) rownames(x@prob)))
		}
	} else {
		children <- NULL
	}	

	glist <- rownames(x@prob)

	for (n in 1:nrow(x@prob)) {
		if (!x@pruned[n] & (x@leaf[n] || (!is.null(clist) & !glist[n] %in% children))) {

			## Calling the gain function
			x@pruned[n] <- !do.call(gain, args=list(p1=x@prob[n,], p2=parent@prob[glist[n],], N=x@n[n], C=C))
		}
	}

	return(x)
}


node.prune <- function(x) { 
	x@pruned[] <- TRUE 
	return(x)
}

node.nmin <- function(x,nmin) {
	x@pruned <- x@n < nmin
	return(x)
}


node.leaf <- function(x) { 
	x@leaf[] <- TRUE 
	return(x)
}

is.leaf <- function(x) { 
	return(x@leaf)
}

is.pruned <- function(x) { 
	return(x@pruned)
}

delete.pruned <- function(x) {
	if (any(x@pruned)) {
		idxnp <- !x@pruned
		x@index <- x@index[idxnp,, drop=FALSE]
		x@prob <- x@prob[idxnp,, drop=FALSE]
		x@counts <- x@counts[idxnp,, drop=FALSE]
		x@pruned <- x@pruned[idxnp,,drop=FALSE]
		x@leaf <- x@leaf[idxnp,,drop=FALSE]
		x@n <- x@n[idxnp,,drop=FALSE]
	}
	return(x)
}


## Used by subtree
select.segment <- function(x, group, position) {
	select <- matrix(FALSE, nrow=nrow(x@index), ncol=ncol(x@index), dimnames=dimnames(x@index))

	## group selection
	if (!is.null(group)) {
		idx.group <- which(x@index[,"group"]==group)
		select[idx.group,"group"] <- TRUE
	} else {
		select[,"group"] <- TRUE
	}

	## position selection
	if (!is.null(position)) {
		idx.pos <- which(x@index[, "position"]==position)
		select[idx.pos,"position"] <- TRUE
	} else {
		select[,"position"] <- TRUE
	}

	## final selection
	idx <- which(rowSums(select)==2)

	if (length(idx)>0) {
		x@index <- x@index[idx,,drop=FALSE]
		x@prob <- x@prob[idx,, drop=FALSE]
		x@counts <- x@counts[idx,, drop=FALSE]
		x@pruned <- x@pruned[idx,, drop=FALSE]
		x@leaf <- x@leaf[idx,, drop=FALSE]
		x@n <- x@n[idx,, drop=FALSE]
	} else {
		x <- NULL
	}

	return(x)
}

## Mining for state prob
node.mine <- function(x, pmin, pmax, state) {
	if (!missing(pmin)) {
		tmp <- rowSums(x@prob[, state, drop=FALSE]) >= pmin
	} else if (!missing(pmax)) {
		tmp <- rowSums(x@prob[, state, drop=FALSE]) < pmax 
	}

	if (sum(tmp)>0) {
		res <- new("cprobd", x@prob[tmp, , drop=FALSE], context=x@path)
	} else {
		res <- NULL
	}	

	return(res)
}


## Merging two nodes
node.merge <- function(x, y, segmented) {
	if (segmented) {
		x@prob <- rbind(x@prob, y@prob)
		x@counts <- rbind(x@counts, y@counts)
		x@n <- c(x@n, y@n)
		x@pruned <- c(x@pruned, y@pruned)
		x@leaf <- c(x@leaf, y@leaf)
	} else {
		x@counts <- rbind(x@prob, y@counts)
		x@prob <- x@counts/sum(x@counts)
		x@n <- x@n+y@n
		x@pruned <- x@pruned & y@pruned
		x@leaf <- x@leaf & y@leaf
	}

	return(x)
}

## set leave
set.leaves <- function(x, clist, cplist) {
	child <- which(x@path==cplist)
	pchild <- unique(as.vector(unlist(clist[child])))
	leaves <- !rownames(x@leaf) %in% pchild
	if (sum(leaves)>0) { x@leaf[leaves] <- TRUE }

	return(x)
}


which.child <- function(PST) {
	class <- class(PST)
	child.list <- names(PST)[unlist(lapply(PST, is, class))]
	return(child.list)
}




