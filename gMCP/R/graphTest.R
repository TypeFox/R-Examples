#' Multiple testing using graphs
#' 
#' Implements the graphical test procedure described in Bretz et al. (2009).
#' Note that the gMCP function in the gMCP package performs the same task.
#' 
#' 
#' @param pvalues Either a vector or a matrix containing the local p-values for
#' the hypotheses in the rows.
#' @param weights Initial weight levels for the test procedure, in case of
#' multiple graphs this needs to be a matrix.
#' @param alpha Overall alpha level of the procedure. For entangled graphs
#' \code{alpha} should be a numeric vector of length equal to the number of
#' graphs, each element specifying the partial alpha for the respective graph.
#' The overall alpha level equals \code{sum(alpha)}.
#' @param G For simple graphs \code{G} should be a numeric matrix determining
#' the graph underlying the test procedure. Note that the diagonal need to
#' contain only 0s, while the rows need to sum to 1.  For entangled graphs it
#' needs to be a list containing the different graph matrices as elements.
#' @param graph As an alternative to the specification via \code{weights} and
#' \code{G} one can also hand over a \code{graphMCP} object to the code.
#' \code{graphMCP} objects can be created for example with the \code{graphGUI}
#' function.
#' @param cr Correlation matrix that should be used for the parametric test.
#' If \code{cr==NULL} the Bonferroni based test procedure is used.
#' @param verbose If verbose is TRUE, additional information about the
#' graphical rejection procedure is displayed.
#' @param test In the parametric case there is more than one way to handle
#' subgraphs with less than the full alpha. If the parameter \code{test} is
#' missing, the tests are performed as described by Bretz et al. (2011), i.e.
#' tests of intersection null hypotheses always exhaust the full alpha level
#' even if the sum of weights is strictly smaller than one. If
#' \code{test="simple-parametric"} the tests are performed as defined in
#' Equation (3) of Bretz et al. (2011).
#' @param upscale Logical. If \code{upscale=FALSE} then for each intersection 
#' of hypotheses (i.e. each subgraph) a weighted test is performed at the 
#' possibly reduced level alpha of sum(w)*alpha, 
#' where sum(w) is the sum of all node weights in this subset.
#' If \code{upscale=TRUE} all weights are upscaled, so that sum(w)=1.
#' @return A vector or a matrix containing the test results for the hypotheses
#' under consideration. Significant tests are denoted by a 1, non-significant
#' results by a 0.
#' @references
#' 
#' Bretz, F., Maurer, W., Brannath, W. and Posch, M. (2009) A graphical
#' approach to sequentially rejective multiple test procedures. Statistics in
#' Medicine, 28, 586--604
#' 
#' Bretz, F., Maurer, W. and Hommel, G. (2010) Test and power considerations
#' for multiple endpoint analyses using sequentially rejective graphical
#' procedures, to appear in Statistics in Medicine
#' @keywords htest
#' @examples
#' 
#' 
#' #### example from Bretz et al. (2010)
#' weights <- c(1/3, 1/3, 1/3, 0, 0, 0)
#' graph <- rbind(c(0,       0.5, 0,     0.5, 0,      0),
#'                c(1/3,     0,   1/3,    0,   1/3,    0),
#'                c(0,       0.5, 0,     0,   0,      0.5),
#'                c(0,       1,   0,     0,   0,      0),
#'                c(0.5,     0,   0.5,   0,   0,      0),
#'                c(0,       1,   0,     0,   0,      0))
#' pvals <- c(0.1, 0.008, 0.005, 0.15, 0.04, 0.006)
#' graphTest(pvals, weights, alpha=0.025, graph)
#' 
#' ## observe graphical procedure in detail
#' graphTest(pvals, weights, alpha=0.025, graph, verbose = TRUE)
#' 
#' ## now use many p-values (useful for power simulations)
#' pvals <- matrix(rbeta(6e4, 1, 30), ncol = 6)
#' out <- graphTest(pvals, weights, alpha=0.025, graph)
#' head(out)
#' 
#' ## example using multiple graphs (instead of 1)
#' G1 <- rbind(c(0,0.5,0.5,0,0), c(0,0,1,0,0),
#'             c(0, 0, 0, 1-0.01, 0.01), c(0, 1, 0, 0, 0),
#'             c(0, 0, 0, 0, 0))
#' G2 <- rbind(c(0,0,1,0,0), c(0.5,0,0.5,0,0),
#'             c(0, 0, 0, 0.01, 1-0.01), c(0, 0, 0, 0, 0),
#'             c(1, 0, 0, 0, 0))
#' weights <- rbind(c(1, 0, 0, 0, 0), c(0, 1, 0, 0, 0))
#' pvals <- c(0.012, 0.025, 0.005, 0.0015, 0.0045)
#' out <- graphTest(pvals, weights, alpha=c(0.0125, 0.0125), G=list(G1, G2), verbose = TRUE)
#' 
#' ## now again with many p-values
#' pvals <- matrix(rbeta(5e4, 1, 30), ncol = 5)
#' out <- graphTest(pvals, weights, alpha=c(0.0125, 0.0125), G=list(G1, G2))
#' head(out)
#' 
#' @export graphTest

graphTest <- function(pvalues, weights = NULL, alpha = 0.05, G = NULL, cr = NULL, graph = NULL, verbose = FALSE, test, upscale=FALSE) {
	usegraph <- !is.null(graph)
	if (!is.list(G) && length(alpha)!=1) {
		stop("Length of 'alpha' should be one for only one graph.")
	} else if (is.list(G) && length(alpha)!=length(G)) {
		stop("length of 'alpha' and 'G' should match")
	}
	if(usegraph & !(class(graph) %in% c("graphMCP", "entangledMCP")))
		stop("Parameter 'graph' needs to be an object of class 'graphMCP' or 'entangledMCP'.")
	if(usegraph & ((!is.null(weights)||!is.null(G))))
		stop("If 'graph' is specified, don't set 'weights' or 'G'.")
	if(usegraph) {    
		## get alpha vector and transition matrix/matrices
    if ("graphMCP" %in% class(graph)) {
		  alphas <- graph@weights * alpha
		  weights <- graph@weights
		  G <- graph@m
    } else {
      G <- getMatrices(graph)
      alpha <-alpha*graph@weights         # Spread alpha between graphs
      alphas <- getWeights(graph) * alpha # Spread alphas between hypotheses
    }
	} else {
		alphas <- weights * alpha
	}
	
	nH <- ifelse(!is.matrix(pvalues), length(pvalues), ncol(pvalues))
	nH <- as.integer(nH)
	checkArgs(pvalues, alphas, G, nH)
	
	if (!is.null(cr)) { # parametric case
		if (is.list(G)) {
			stop("The parametric case does not support multiple graphs yet.")
		}
		hint <- generateWeights(G, weights)
		out <- matrix(0, nrow=0, ncol=dim(pvalues)[2])
		colnames(out) <- colnames(G)
		for (i in 1:(dim(pvalues)[1])) {
			adjP <- generatePvals(G, weights, cr, pvalues[i,], hint=hint, upscale=upscale)
			out <- rbind(out, ifelse(adjP<=alpha,1,0))
		}
		return(out)
	} else { # non-parametric case		
		if(is.list(G) || upscale==FALSE){
		  if(is.list(G)) {
		    nGraphs <- length(G)
		    G <- c(unlist(G))
		  } else {
		    nGraphs <- 1
		  }
			if (upscale==TRUE) stop("Upscale=TRUE with list of graphs currently not supported")
			if(!is.matrix(pvalues)){
			  res <- .C("graphproc", h=double(nH), a=as.double(alphas), G=as.double(G),
			            as.double(pvalues), nH, as.double(G), as.integer(nGraphs),
			            as.integer(verbose), as.integer(upscale))
			  out <- c(H = res$h)
			  attr(out, "last.alphas") <- res$a
			  attr(out, "last.G") <- matrix(res$G, ncol = nH)
			  return(out)
			} else {
			  nCount <- as.integer(nrow(pvalues))
			  res <- .C("graphmult", h=double(nH*nCount), double(nH),
			            as.double(alphas), double(nGraphs*nH),
			            as.double(G), as.double(G), as.double(G),
			            as.double(pvalues), double(nH), nCount, nH,
			            as.integer(nGraphs), as.integer(verbose))
			  out <- matrix(res$h, nrow = nCount)
			  if(is.null(colnames(G))) {
			    colnames(out) <- paste("H", 1:nH, sep="")
			  } else {
			    colnames(out) <- colnames(G)
			  }
			  return(out)
			}
		}
	  out <- matrix(0, nrow=0, ncol=dim(pvalues)[2])
	  colnames(out) <- colnames(G)
	  for (i in 1:(dim(pvalues)[1])) {
	    adjP <- gMCP(new("graphMCP", m=G, weights=weights), pvalues[i,], upscale=upscale)@adjPValues
	    out <- rbind(out, ifelse(adjP<=alpha,1,0))
	  }
	  return(out)
	}
}

checkArgs <- function(pvalues, alphas, G, nH){
	# debug print(pvalues); print(alphas); print(G); print(nH)
	alplen <- ifelse(is.matrix(alphas), ncol(alphas),
			length(alphas))
	if(any(alphas < 0) | any(alphas > 1))
		stop("entries of alphas need to be in [0,1].")
	if(any(pvalues < 0) | any(pvalues > 1))
		stop("entries of pvalues need to be in [0,1].")
	if(is.list(G)){
		nGraphs <- length(G)
		for(i in 1:nGraphs){
			if(any(G[[i]] < 0) | any(G[[i]] > 1))
				stop("entries of G need to be in [0,1]")
			if(any(rowSums(G[[i]]) > 1))
				stop("rows of G need to sum to values <= 1")
			if(nrow(G[[i]]) != ncol(G[[i]]))
				stop("non-quadratic matrix G")
			if(nH != nrow(G[[i]]))
				stop("non-conforming pvalues and G.")
			if(alplen != nrow(G[[i]]))
				stop("non-conforming alphas and G.")
		}
		if(!is.matrix(alphas)){
			stop("alphas needs to have as many rows as there are graphs")
		} else {
			if(nrow(alphas) != nGraphs)
				stop("alphas needs to have as many rows as there are graphs")
		}
		
	} else {
    if(!is.matrix(G) || !is.numeric(G)) {
      stop("G must be a numeric matrix")
    }
		if(any(G < 0) | any(G > 1))
			stop("entries of G need to be in [0,1]")
		if(any(rowSums(G) > 1))
			stop("rows of G need to sum to values <= 1")
		if(nrow(G) != ncol(G))
			stop("non-quadratic matrix G") 
		if(nH != nrow(G))
			stop("non-conforming pvalues and G.")
		if(alplen != nrow(G))
			stop("non-conforming alphas and G.")
	}
}

convert <- function(g){
	## converts a graph object as
	## used in the gMCP library
	## and returns vector of alphas
	## and significance levels
	if(class(g) != "graphMCP")
		stop("g needs to an object of class graphMCP")
	Hnams <- g@nodes
	nH <- length(Hnams)
	alphas <- numeric(nH)
	G <- matrix(0, nrow = nH, ncol = nH)
	for(i in 1:nH){
		alphas[i] <- g@nodeAttr@data[[i]]$alpha
	}
	names(alphas) <- Hnams
	
	nams <- names(g@edgeAttr)
	for(nam in nams){
		nam2 <- strsplit(nam, "\\|")
		indx <- grep(nam2[[1]][1], Hnams)
		indy <- grep(nam2[[1]][2], Hnams)
		wgt <- g@edgeAttr@data[[nam]]$weight
		G[indx, indy] <- wgt
	}
	dimnames(G) <- list(Hnams, Hnams)
	list(alphas=alphas, G=G)
}
