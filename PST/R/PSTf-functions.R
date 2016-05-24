likelihood <- function(object, log=TRUE) {
	debut.lik <- Sys.time()
	message(" [>] computing sequence(s) likelihood ...", appendLF=FALSE)
	lik <- suppressMessages(predict(object, object@data, object@cdata, group=object@group, decomp=TRUE))
	if (log) { lik <- sum(log(lik)) }
	fin.lik <- Sys.time()
	message(" (", format(round(fin.lik-debut.lik, 3)), ")")
	return(lik)
} 



node.list <- function(x, pruned=FALSE) {

	plist <- lapply(x, names)
	plist
}

pruned.nodes <- function(x) {
	pruned.list <- lapply(x, is.pruned)
	pruned.list <- unlist(pruned.list)
	return(pruned.list)
}

leaves.list <- function(x) {
	res <- lapply(x, function(x) { all(x@leaf) } )
	res <- unlist(res)
	return(res)
}


leaves.count <- function(x) {
	res <- lapply(x, function(x) { sum(x@leaf) } )
	res <- unlist(res)
	return(res)
}


nodes.count <- function(x) {
	res <- lapply(x, function(x) { sum(!x@leaf) } )
	res <- unlist(res)
	return(res)
}

## list of leaves
leaves <- function(x) {
	res <- NULL
	for (i in 1:length(x)) {
		nodes.names <- node.list(x[[i]]) 
		tmp <- lapply(x[[i]], function(x) { x@leaf } )
		tmp <- unlist(tmp)
		res <- c(res, nodes.names[tmp])
	}
	return(res)
}


is.stationary <- function(x) {
	all(is.na(x[[1]][["e"]]@index[, "position"]))
}

has.cdata <- function(x) {
	(nrow(x@cdata) > 0)
}


