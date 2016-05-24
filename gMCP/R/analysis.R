#' Analysis of a gMCP-Graph
#' 
#' Creates LaTeX code that represents the given graph.
#' 
#' In the moment it is only tested whether each node is accessible from each
#' other node. Further analysis will be added in future versions.
#' 
#' @param graph A graph of class \code{\link{graphMCP}}.
#' @param file A connection, or a character string naming the file to print to.
#' If \code{""} (the default), the analysis is printed to the standard output
#' connection (the console unless redirected by \code{sink}).  If it is
#' \code{"|cmd"}, the output is piped to the command given by \code{cmd}, by
#' opening a pipe connection [taken from the manual page of \code{cat}, which
#' is called in this function].
#' @return A character string that contains the printed analysis.
#' @author Kornelius Rohmeyer \email{rohmeyer@@small-projects.de}
#' @seealso \code{\link{graphMCP}}
#' @keywords graphs
#' @examples
#' 
#' 
#' g <- BonferroniHolm(5)
#' 
#' graphAnalysis(g)
#' 
#' 
#' @export graphAnalysis
graphAnalysis <- function(graph, file="") {
	if ("entangledMCP" %in% class(graph)) {
		return("Analysis of entangled graphs is not yet supported.")
	}
	result <- checkOptimal(graph, verbose=FALSE)
	cat(result, file=file)
	return(invisible(result))
}

accessible <- function(graph, node) {
	m <- graph@m
	#m <- ifelse(!is.na(as.num(m)) && as.num(m) != 0, 0, 1)
	s <- which(getNodes(graph)==node)
	n <- dim(m)[1]
	ac <- rep(FALSE, n)
	followed <- rep(FALSE, n)
	for (j in 1:n) {
		for (i in 1:n) {
			if (followed[i]==FALSE && (ac[i]==TRUE || i == s)) {
				ac[m[i,]!=0] <- TRUE
				followed[i] <- TRUE
			}
		}
	}
	return(ac)
}

# Checks whether the weight of all subgraphs sums up to 1.
# 
# Checks whether the weight of all subgraphs sums up to 1.
# 
# @param g Graph of class \code{\link{graphMCP}} (TODO: or \code{\link{entangledMCP}}).
# @param verbose Logical whether to print the return value.
# @author Kornelius Rohmeyer \email{rohmeyer@@small-projects.de}
# @keywords Character string containing the analysis.
# @examples
#
# graph <- HungEtWang2010()
# # Does not work with variables:
# gMCP:::checkOptimal(graph)
# graph <- replaceVariables(graph, list("tau"=0.5,"omega"=0.5, "nu"=0.5))
# gMCP:::checkOptimal(graph)
#
checkOptimal <- function (graph, verbose=TRUE) {  
  if (any(is.na(parse2numeric(graph, force=TRUE)@m))) return("Graph contains variables (or unparsable terms). No analysis possible.")
	nodes <- getNodes(graph)[getWeights(graph)!=0]
	s <- ""
	n <- dim(graph@m)[1]
	all.weights <- generateWeights(graph@m,getWeights(graph))[,c((n+1):(2*n))]
	if (!isTRUE(all.equal(rowSums(all.weights), rep(1, dim(all.weights)[1])))) {
	    s <- "The graph is not optimal."
    if (isTRUE(all.equal(sum(getWeights(graph)),1))) {
      wrong <- unlist(lapply(rowSums(all.weights), function(x) { !isTRUE(all.equal(x, 1)) }))
      s <- paste(s, "\nThe weights in ",sum(wrong)," subgraphs do not sum up to 1.", sep="")
    } else {
      s <- paste(s, "The weights of the original graph do not sum up to 1.", sep="\n")
    }
	} else {
	  s <- "For all subgraphs the weights sum up to 1."
	}
	#for (n in nodes) {
	#	notAccessible <- setdiff(getNodes(graph)[!accessible(graph, n)], n)
	#	if (length(notAccessible)>0) {
	#		if (s=="") {
	#			s <- "The graph is not optimal.\nBy adding edges the test can be improved uniformly."
	#			#s <- paste(s, "Or set exhaustAlpha=TRUE to do an alpha exhaustive test as described in Bretz et al. (2011).\n", sep="\n")
	#		}
	#		s <- paste(s, "There is no path from node ",n, " to ", paste(notAccessible, collapse=", "), "\n", sep="")
	#	}
	#}
	#if (s=="") s <- "From each node with positive weight paths exist to all other nodes.\n"
	if (verbose) cat(s)
	return(s)
}

as.num <- function(x) {
	options(warn=-1)
	result <- as.numeric(x)
	options(warn=0)
	return(result)
}