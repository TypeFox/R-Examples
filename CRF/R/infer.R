#' Inference method for small graphs
#' 
#' Computing the partition function and marginal probabilities
#' 
#' Exact inference for small graphs with brute-force counting 
#' 
#' @param crf The CRF
#' @return This function will return a list with components:
#'   \item{node.bel}{Node belief. It is a matrix with \code{crf$n.nodes} rows and \code{crf$max.state} columns.}
#'   \item{edge.bel}{Edge belief. It is a list of matrices. The size of list is \code{crf$n.edges} and 
#'     the matrix \code{i} has \code{crf$n.states[crf$edges[i,1]]} rows and \code{crf$n.states[crf$edges[i,2]]} columns.}
#'   \item{logZ}{The logarithmic value of CRF normalization factor Z.}
#' 
#' @examples
#' 
#' library(CRF)
#' data(Small)
#' i <- infer.exact(Small$crf)
#' 
#' @export
infer.exact <- function(crf)
	.Call(Infer_Exact, crf)



#' Inference method for chain-structured graphs
#' 
#' Computing the partition function and marginal probabilities
#' 
#' Exact inference for chain-structured graphs with the forward-backward algorithm 
#' 
#' @param crf The CRF
#' @return This function will return a list with components:
#'   \item{node.bel}{Node belief. It is a matrix with \code{crf$n.nodes} rows and \code{crf$max.state} columns.}
#'   \item{edge.bel}{Edge belief. It is a list of matrices. The size of list is \code{crf$n.edges} and 
#'     the matrix \code{i} has \code{crf$n.states[crf$edges[i,1]]} rows and \code{crf$n.states[crf$edges[i,2]]} columns.}
#'   \item{logZ}{The logarithmic value of CRF normalization factor Z.}
#' 
#' @examples
#' 
#' library(CRF)
#' data(Small)
#' i <- infer.chain(Small$crf)
#' 
#' @export
infer.chain <- function(crf)
	.Call(Infer_Chain, crf)



#' Inference method for tree- and forest-structured graphs
#' 
#' Computing the partition function and marginal probabilities
#' 
#' Exact inference for tree- and forest-structured graphs with sum-product belief propagation 
#' 
#' @param crf The CRF
#' @return This function will return a list with components:
#'   \item{node.bel}{Node belief. It is a matrix with \code{crf$n.nodes} rows and \code{crf$max.state} columns.}
#'   \item{edge.bel}{Edge belief. It is a list of matrices. The size of list is \code{crf$n.edges} and 
#'     the matrix \code{i} has \code{crf$n.states[crf$edges[i,1]]} rows and \code{crf$n.states[crf$edges[i,2]]} columns.}
#'   \item{logZ}{The logarithmic value of CRF normalization factor Z.}
#' 
#' @examples
#' 
#' library(CRF)
#' data(Small)
#' i <- infer.tree(Small$crf)
#' 
#' @export
infer.tree <- function(crf)
	.Call(Infer_Tree, crf)



#' Conditional inference method
#' 
#' Computing the partition function and marginal probabilities
#' 
#' Conditional inference (takes another inference method as input) 
#' 
#' @param crf The CRF
#' @param clamped The vector of fixed values for clamped nodes, 0 for unfixed nodes
#' @param infer.method The inference method to solve the clamped CRF
#' @param ... The parameters for \code{infer.method}
#' @return This function will return a list with components:
#'   \item{node.bel}{Node belief. It is a matrix with \code{crf$n.nodes} rows and \code{crf$max.state} columns.}
#'   \item{edge.bel}{Edge belief. It is a list of matrices. The size of list is \code{crf$n.edges} and 
#'     the matrix \code{i} has \code{crf$n.states[crf$edges[i,1]]} rows and \code{crf$n.states[crf$edges[i,2]]} columns.}
#'   \item{logZ}{The logarithmic value of CRF normalization factor Z.}
#' 
#' @examples
#' 
#' library(CRF)
#' data(Small)
#' i <- infer.conditional(Small$crf, c(0,1,0,0), infer.exact)
#' 
#' @export
infer.conditional <- function(crf, clamped, infer.method, ...)
{
	belief <- list()
	belief$node.bel <- array(0, dim=c(crf$n.nodes, crf$max.state))
	belief$edge.bel <- lapply(1:crf$n.edges, function(i) array(0, dim=c(crf$n.states[crf$edges[i,1]], crf$n.states[crf$edges[i,2]])))
	newcrf <- clamp.crf(crf, clamped)
	b <- infer.method(newcrf, ...)
	belief$node.bel[newcrf$node.id, 1:newcrf$max.state] <- b$node.bel
	belief$edge.bel[newcrf$edge.id] <- b$edge.bel
	belief$logZ <- b$logZ
	belief$node.bel[cbind(which(clamped != 0), clamped[clamped != 0])] <- 1
	e <- newcrf$original$edges
	e0 <- which(clamped[e[,1]] != 0 & clamped[e[,2]] != 0)
	e1 <- which(clamped[e[,1]] != 0 & clamped[e[,2]] == 0)
	e2 <- which(clamped[e[,1]] == 0 & clamped[e[,2]] != 0)
	for (i in e0) belief$edge.bel[[i]][clamped[e[i,1]], clamped[e[i,2]]] <- 1
	for (i in e1) belief$edge.bel[[i]][clamped[e[i,1]],] <- b$node.bel[newcrf$node.map[e[i,2]], 1:crf$n.states[e[i,2]]]
	for (i in e2) belief$edge.bel[[i]][,clamped[e[i,2]]] <- b$node.bel[newcrf$node.map[e[i,1]], 1:crf$n.states[e[i,1]]]
	belief
}



#' Inference method for graphs with a small cutset
#' 
#' Computing the partition function and marginal probabilities
#' 
#' Exact inference for graphs with a small cutset using cutset conditioning 
#' 
#' @param crf The CRF
#' @param cutset A vector of nodes in the cutset
#' @param engine The underlying engine for cutset decoding, possible values are "default", "none", "exact", "chain", and "tree".
#' @return This function will return a list with components:
#'   \item{node.bel}{Node belief. It is a matrix with \code{crf$n.nodes} rows and \code{crf$max.state} columns.}
#'   \item{edge.bel}{Edge belief. It is a list of matrices. The size of list is \code{crf$n.edges} and 
#'     the matrix \code{i} has \code{crf$n.states[crf$edges[i,1]]} rows and \code{crf$n.states[crf$edges[i,2]]} columns.}
#'   \item{logZ}{The logarithmic value of CRF normalization factor Z.}
#' 
#' @examples
#' 
#' library(CRF)
#' data(Small)
#' i <- infer.cutset(Small$crf, c(2))
#' 
#' @export
infer.cutset <- function(crf, cutset, engine = "default")
{
	engine.id <- c("default"=-1, "none"=0, "exact"=1, "chain"=2, "tree"=3);
	clamped <- rep(0, crf$n.nodes)
	clamped[cutset] <- 1
	newcrf <- clamp.crf(crf, clamped)
	.Call(Infer_Cutset, newcrf, engine.id[engine])
}



#' Inference method for low-treewidth graphs
#' 
#' Computing the partition function and marginal probabilities
#' 
#' Exact decoding for low-treewidth graphs using junction trees 
#' 
#' @param crf The CRF
#' @return This function will return a list with components:
#'   \item{node.bel}{Node belief. It is a matrix with \code{crf$n.nodes} rows and \code{crf$max.state} columns.}
#'   \item{edge.bel}{Edge belief. It is a list of matrices. The size of list is \code{crf$n.edges} and 
#'     the matrix \code{i} has \code{crf$n.states[crf$edges[i,1]]} rows and \code{crf$n.states[crf$edges[i,2]]} columns.}
#'   \item{logZ}{The logarithmic value of CRF normalization factor Z.}
#' 
#' @examples
#' 
#' library(CRF)
#' data(Small)
#' i <- infer.junction(Small$crf)
#' 
#' @export
infer.junction <- function(crf)
	.Call(Infer_Junction, crf)



#' Inference method using sampling
#' 
#' Computing the partition function and marginal probabilities
#' 
#' Approximate inference using sampling (takes a sampling method as input) 
#' 
#' @param crf The CRF
#' @param sample.method The sampling method
#' @param ... The parameters for \code{sample.method}
#' @return This function will return a list with components:
#'   \item{node.bel}{Node belief. It is a matrix with \code{crf$n.nodes} rows and \code{crf$max.state} columns.}
#'   \item{edge.bel}{Edge belief. It is a list of matrices. The size of list is \code{crf$n.edges} and 
#'     the matrix \code{i} has \code{crf$n.states[crf$edges[i,1]]} rows and \code{crf$n.states[crf$edges[i,2]]} columns.}
#'   \item{logZ}{The logarithmic value of CRF normalization factor Z.}
#' 
#' @examples
#' 
#' library(CRF)
#' data(Small)
#' i <- infer.sample(Small$crf, sample.exact, 10000)
#' 
#' @export
infer.sample <- function(crf, sample.method, ...)
	.Call(Infer_Sample, crf, sample.method(crf, ...))



#' Inference method using loopy belief propagation
#' 
#' Computing the partition function and marginal probabilities
#' 
#' Approximate inference using sum-product loopy belief propagation 
#' 
#' @param crf The CRF
#' @param max.iter The maximum allowed iterations of termination criteria
#' @param cutoff The convergence cutoff of termination criteria
#' @param verbose Non-negative integer to control the tracing informtion in algorithm
#' @return This function will return a list with components:
#'   \item{node.bel}{Node belief. It is a matrix with \code{crf$n.nodes} rows and \code{crf$max.state} columns.}
#'   \item{edge.bel}{Edge belief. It is a list of matrices. The size of list is \code{crf$n.edges} and 
#'     the matrix \code{i} has \code{crf$n.states[crf$edges[i,1]]} rows and \code{crf$n.states[crf$edges[i,2]]} columns.}
#'   \item{logZ}{The logarithmic value of CRF normalization factor Z.}
#' 
#' @examples
#' 
#' library(CRF)
#' data(Small)
#' i <- infer.lbp(Small$crf)
#' 
#' @export
infer.lbp <- function(crf, max.iter = 10000, cutoff = 1e-4, verbose = 0)
	.Call(Infer_LBP, crf, max.iter, cutoff, verbose)



#' Inference method using tree-reweighted belief propagation
#' 
#' Computing the partition function and marginal probabilities
#' 
#' Approximate inference using sum-product tree-reweighted belief propagation
#' 
#' @param crf The CRF
#' @param max.iter The maximum allowed iterations of termination criteria
#' @param cutoff The convergence cutoff of termination criteria
#' @param verbose Non-negative integer to control the tracing informtion in algorithm
#' @return This function will return a list with components:
#'   \item{node.bel}{Node belief. It is a matrix with \code{crf$n.nodes} rows and \code{crf$max.state} columns.}
#'   \item{edge.bel}{Edge belief. It is a list of matrices. The size of list is \code{crf$n.edges} and 
#'     the matrix \code{i} has \code{crf$n.states[crf$edges[i,1]]} rows and \code{crf$n.states[crf$edges[i,2]]} columns.}
#'   \item{logZ}{The logarithmic value of CRF normalization factor Z.}
#' 
#' @examples
#' 
#' library(CRF)
#' data(Small)
#' i <- infer.trbp(Small$crf)
#' 
#' @export
infer.trbp <- function(crf, max.iter = 10000, cutoff = 1e-4, verbose = 0)
	.Call(Infer_TRBP, crf, max.iter, cutoff, verbose)
