#' Convenience for Clustering
#'
#'	Provides a convenience interface for using MPD-style hostfiles to generate cluster objetcs.
#'	The hostfile should be written as a text file using the MPD style:
#'	one line for each host, which can be followed by a colon and a number
#'	indicating the number of processes to be started on that host. An example
#'	hostfile for starting three processes on two hosts named avalon and
#'	glastonbury would be:
#'
#'	avalon
#'	glastonbury:2
#'
#' @param name Filename of the hostfile.
#' @examples
#'	\dontrun{
#'		library(parallel)
#'		cl = makePSOCKcluster(machinefile("mpd.hosts"))
#'		stopCluster(cl)
#'	}
#' @export
machinefile <- function(name) {
	x <- utils::read.table(name, sep=":", header=FALSE, stringsAsFactors=FALSE, fill=TRUE)
	ret <- c()
	for (i in 1:(dim(x)[1])) { 
		if (is.null(x[i,2]) | is.na (x[i,2])) x[i,2] <- 1
		ret <-c(ret, rep(x[i,1], as.numeric(x[i,2])))
	}
	return(ret)
}

# Internal use only
clusterRun <- function (cl, model, L) {
	N <- dim(L)[1]
	sp <- parallel::clusterSplit(cl, 1:N)
	tmp.res <- parallel::clusterApply(cl, sp, 
			 fun = function(idx, x) model(x[idx,]), L)
	n.outs <- length(unlist(tmp.res))/N
	res <- array( dim=c(N, n.outs));
	for (i in 1:(length(sp))) 
		res[sp[[i]], 1:n.outs] <- t(tmp.res[[i]])
	return (res)
}
