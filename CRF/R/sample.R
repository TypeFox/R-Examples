#' Sampling method for small graphs
#' 
#' Generating samples from the distribution
#' 
#' Exact sampling for small graphs with brute-force inverse cumulative distribution 
#' 
#' @param crf The CRF
#' @param size The sample size
#' @return This function will return a matrix with \code{size} rows and \code{crf$n.nodes} columns,
#'   in which each row is a sampled configuration.
#' 
#' @examples
#' 
#' library(CRF)
#' data(Small)
#' s <- sample.exact(Small$crf, 100)
#' 
#' @export
sample.exact <- function(crf, size)
	.Call(Sample_Exact, crf, size)



#' Sampling method for chain-structured graphs
#' 
#' Generating samples from the distribution
#' 
#' Exact sampling for chain-structured graphs with the forward-filter backward-sample algorithm 
#' 
#' @param crf The CRF
#' @param size The sample size
#' @return This function will return a matrix with \code{size} rows and \code{crf$n.nodes} columns,
#'   in which each row is a sampled configuration.
#' 
#' @examples
#' 
#' library(CRF)
#' data(Small)
#' s <- sample.chain(Small$crf, 100)
#' 
#' @export
sample.chain <- function(crf, size)
	.Call(Sample_Chain, crf, size)



#' Sampling method for tree- and forest-structured graphs
#' 
#' Generating samples from the distribution
#' 
#' Exact sampling for tree- and forest-structured graphs with sum-product belief propagation and backward-sampling 
#' 
#' @param crf The CRF
#' @param size The sample size
#' @return This function will return a matrix with \code{size} rows and \code{crf$n.nodes} columns,
#'   in which each row is a sampled configuration.
#' 
#' @examples
#' 
#' library(CRF)
#' data(Small)
#' s <- sample.tree(Small$crf, 100)
#' 
#' @export
sample.tree <- function(crf, size)
	.Call(Sample_Tree, crf, size)



#' Conditional sampling method
#' 
#' Generating samples from the distribution
#' 
#' Conditional sampling (takes another sampling method as input) 
#' 
#' @param crf The CRF
#' @param size The sample size
#' @param clamped The vector of fixed values for clamped nodes, 0 for unfixed nodes
#' @param sample.method The sampling method to solve the clamped CRF
#' @param ... The parameters for \code{sample.method}
#' @return This function will return a matrix with \code{size} rows and \code{crf$n.nodes} columns,
#'   in which each row is a sampled configuration.
#' 
#' @examples
#' 
#' library(CRF)
#' data(Small)
#' s <- sample.conditional(Small$crf, 100, c(0,1,0,0), sample.exact)
#' 
#' @export
sample.conditional <- function(crf, size, clamped, sample.method, ...)
{
	newcrf <- clamp.crf(crf, clamped)
	s <- sample.method(newcrf, size, ...)
	samples <- matrix(rep(clamped, nrow(s)), nrow=nrow(s), ncol=length(clamped), byrow=TRUE)
	samples[,newcrf$node.id] <- s
	samples
}



#' Sampling method for graphs with a small cutset
#' 
#' Generating samples from the distribution
#' 
#' Exact sampling for graphs with a small cutset using cutset conditioning 
#' 
#' @param crf The CRF
#' @param size The sample size
#' @param cutset A vector of nodes in the cutset
#' @param engine The underlying engine for cutset sampling, possible values are "default", "none", "exact", "chain", and "tree".
#' @return This function will return a matrix with \code{size} rows and \code{crf$n.nodes} columns,
#'   in which each row is a sampled configuration.
#' 
#' @examples
#' 
#' library(CRF)
#' data(Small)
#' s <- sample.cutset(Small$crf, 100, c(2))
#' 
#' @export
sample.cutset <- function(crf, size, cutset, engine = "default")
{
	engine.id <- c("default"=-1, "none"=0, "exact"=1, "chain"=2, "tree"=3);
	clamped <- rep(0, crf$n.nodes)
	clamped[cutset] <- 1
	newcrf <- clamp.crf(crf, clamped)
	.Call(Sample_Cutset, newcrf, size, engine.id[engine])
}



#' Sampling method for low-treewidth graphs
#' 
#' Generating samples from the distribution
#' 
#' Exact sampling for low-treewidth graphs using junction trees 
#' 
#' @param crf The CRF
#' @param size The sample size
#' @return This function will return a matrix with \code{size} rows and \code{crf$n.nodes} columns,
#'   in which each row is a sampled configuration.
#' 
#' @examples
#' 
#' library(CRF)
#' data(Small)
#' s <- sample.junction(Small$crf, 100)
#' 
#' @export
sample.junction <- function(crf, size)
	.Call(Sample_Junction, crf, size)



#' Sampling method using single-site Gibbs sampler
#' 
#' Generating samples from the distribution
#' 
#' Approximate sampling using a single-site Gibbs sampler
#' 
#' @param crf The CRF
#' @param size The sample size
#' @param burn.in The number of samples at the beginning that will be discarded
#' @param start An initial configuration
#' @return This function will return a matrix with \code{size} rows and \code{crf$n.nodes} columns,
#'   in which each row is a sampled configuration.
#' 
#' @examples
#' 
#' library(CRF)
#' data(Small)
#' s <- sample.gibbs(Small$crf, 100)
#' 
#' @export
sample.gibbs <- function(crf, size, burn.in = 1000, start = apply(crf$node.pot, 1, which.max))
	.Call(Sample_Gibbs, crf, size, burn.in, start)
