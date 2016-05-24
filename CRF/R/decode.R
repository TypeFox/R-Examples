#' Decoding method for small graphs
#' 
#' Computing the most likely configuration for CRF
#' 
#' Exact decoding for small graphs with brute-force search 
#' 
#' @param crf The CRF
#' @return This function will return the most likely configuration, which is a vector of length \code{crf$n.nodes}.
#' 
#' @examples
#' 
#' library(CRF)
#' data(Small)
#' d <- decode.exact(Small$crf)
#'
#' @export 
decode.exact <- function(crf)
	.Call(Decode_Exact, crf)



#' Decoding method for chain-structured graphs
#' 
#' Computing the most likely configuration for CRF
#' 
#' Exact decoding for chain-structured graphs with the Viterbi algorithm.
#' 
#' @param crf The CRF
#' @return This function will return the most likely configuration, which is a vector of length \code{crf$n.nodes}.
#' 
#' @examples
#' 
#' library(CRF)
#' data(Small)
#' d <- decode.chain(Small$crf)
#'
#' @export 
decode.chain <- function(crf)
	.Call(Decode_Chain, crf)



#' Decoding method for tree- and forest-structured graphs
#' 
#' Computing the most likely configuration for CRF
#' 
#' Exact decoding for tree- and forest-structured graphs with max-product belief propagation 
#' 
#' @param crf The CRF
#' @return This function will return the most likely configuration, which is a vector of length \code{crf$n.nodes}.
#' 
#' @examples
#' 
#' library(CRF)
#' data(Small)
#' d <- decode.tree(Small$crf)
#'
#' @export 
decode.tree <- function(crf)
	.Call(Decode_Tree, crf)



#' Conditional decoding method
#' 
#' Computing the most likely configuration for CRF
#' 
#' Conditional decoding (takes another decoding method as input) 
#' 
#' @param crf The CRF
#' @param clamped The vector of fixed values for clamped nodes, 0 for unfixed nodes
#' @param decode.method The decoding method to solve clamped CRF
#' @param ... The parameters for \code{decode.method}
#' @return This function will return the most likely configuration, which is a vector of length \code{crf$n.nodes}.
#' 
#' @examples
#' 
#' library(CRF)
#' data(Small)
#' d <- decode.conditional(Small$crf, c(0,1,0,0), decode.exact)
#'
#' @export 
decode.conditional <- function(crf, clamped, decode.method, ...)
{
	newcrf <- clamp.crf(crf, clamped)
	decode <- clamped
	decode[newcrf$node.id] <- decode.method(newcrf, ...)
	decode
}



#' Decoding method for graphs with a small cutset
#' 
#' Computing the most likely configuration for CRF
#' 
#' Exact decoding for graphs with a small cutset using cutset conditioning 
#' 
#' @param crf The CRF
#' @param cutset A vector of nodes in the cutset
#' @param engine The underlying engine for cutset decoding, possible values are "default", "none", "exact", "chain", and "tree".
#' @param start An initial configuration, a good start will significantly reduce the seraching time
#' @return This function will return the most likely configuration, which is a vector of length \code{crf$n.nodes}.
#' 
#' @examples
#' 
#' library(CRF)
#' data(Small)
#' d <- decode.cutset(Small$crf, c(2))
#'
#' @export 
decode.cutset <- function(crf, cutset, engine = "default", start = apply(crf$node.pot, 1, which.max))
{
	engine.id <- c("default"=-1, "none"=0, "exact"=1, "chain"=2, "tree"=3);
	clamped <- rep(0, crf$n.nodes)
	clamped[cutset] <- 1
	newcrf <- clamp.crf(crf, clamped)
	.Call(Decode_Cutset, newcrf, engine.id[engine], start)
}



#' Decoding method for low-treewidth graphs
#' 
#' Computing the most likely configuration for CRF
#' 
#' Exact decoding for low-treewidth graphs using junction trees 
#' 
#' @param crf The CRF
#' @return This function will return the most likely configuration, which is a vector of length \code{crf$n.nodes}.
#' 
#' @examples
#' 
#' library(CRF)
#' data(Small)
#' d <- decode.junction(Small$crf)
#'
#' @export 
decode.junction <- function(crf)
	.Call(Decode_Junction, crf);



#' Decoding method using sampling
#' 
#' Computing the most likely configuration for CRF
#' 
#' Approximate decoding using sampling (takes a sampling method as input) 
#' 
#' @param crf The CRF
#' @param sample.method The sampling method
#' @param ... The parameters for \code{sample.method}
#' @return This function will return the most likely configuration, which is a vector of length \code{crf$n.nodes}.
#' 
#' @examples
#' 
#' library(CRF)
#' data(Small)
#' d <- decode.sample(Small$crf, sample.exact, 10000)
#'
#' @export 
decode.sample <- function(crf, sample.method, ...)
	.Call(Decode_Sample, crf, sample.method(crf, ...))



#' Decoding method using inference
#' 
#' Computing the most likely configuration for CRF
#' 
#' Approximate decoding using inference (takes an inference method as input) 
#' 
#' @param crf The CRF
#' @param infer.method The inference method
#' @param ... The parameters for \code{infer.method}
#' @return This function will return the most likely configuration, which is a vector of length \code{crf$n.nodes}.
#' 
#' @examples
#' 
#' library(CRF)
#' data(Small)
#' d <- decode.marginal(Small$crf, infer.exact)
#'
#' @export 
decode.marginal <- function(crf, infer.method, ...)
	apply(infer.method(crf, ...)$node.bel, 1, which.max)



#' Decoding method using loopy belief propagation
#' 
#' Computing the most likely configuration for CRF
#' 
#' Approximate decoding using max-product loopy belief propagation 
#' 
#' @param crf The CRF
#' @param max.iter The maximum allowed iterations of termination criteria
#' @param cutoff The convergence cutoff of termination criteria
#' @param verbose Non-negative integer to control the tracing informtion in algorithm
#' @return This function will return the most likely configuration, which is a vector of length \code{crf$n.nodes}.
#' 
#' @examples
#' 
#' library(CRF)
#' data(Small)
#' d <- decode.lbp(Small$crf)
#'
#' @export 
decode.lbp <- function(crf, max.iter = 10000, cutoff = 1e-4, verbose = 0)
	.Call(Decode_LBP, crf, max.iter, cutoff, verbose)



#' Decoding method using tree-reweighted belief propagation
#' 
#' Computing the most likely configuration for CRF
#' 
#' Approximate decoding using max-product tree-reweighted belief propagtion 
#' 
#' @param crf The CRF
#' @param max.iter The maximum allowed iterations of termination criteria
#' @param cutoff The convergence cutoff of termination criteria
#' @param verbose Non-negative integer to control the tracing informtion in algorithm
#' @return This function will return the most likely configuration, which is a vector of length \code{crf$n.nodes}.
#' 
#' @examples
#' 
#' library(CRF)
#' data(Small)
#' d <- decode.trbp(Small$crf)
#'
#' @export 
decode.trbp <- function(crf, max.iter = 10000, cutoff = 1e-4, verbose = 0)
	.Call(Decode_TRBP, crf, max.iter, cutoff, verbose)



#' Decoding method using greedy algorithm
#' 
#' Computing the most likely configuration for CRF
#' 
#' Approximate decoding with greedy algorithm 
#' 
#' @param crf The CRF
#' @param restart Non-negative integer to control how many restart iterations are repeated
#' @param start An initial configuration, a good start will significantly reduce the seraching time
#' @return This function will return the most likely configuration, which is a vector of length \code{crf$n.nodes}.
#' 
#' @examples
#' 
#' library(CRF)
#' data(Small)
#' d <- decode.greedy(Small$crf)
#'
#' @export 
decode.greedy <- function(crf, restart = 0, start = apply(crf$node.pot, 1, which.max))
	.Call(Decode_Greedy, crf, restart, start)



#' Decoding method using iterated conditional modes algorithm
#' 
#' Computing the most likely configuration for CRF
#' 
#' Approximate decoding with the iterated conditional modes algorithm 
#' 
#' @param crf The CRF
#' @param restart Non-negative integer to control how many restart iterations are repeated
#' @param start An initial configuration, a good start will significantly reduce the seraching time
#' @return This function will return the most likely configuration, which is a vector of length \code{crf$n.nodes}.
#' 
#' @examples
#' 
#' library(CRF)
#' data(Small)
#' d <- decode.icm(Small$crf)
#'
#' @export 
decode.icm <- function(crf, restart = 0, start = apply(crf$node.pot, 1, which.max))
	.Call(Decode_ICM, crf, restart, start)



#' Decoding method using block iterated conditional modes algorithm
#' 
#' Computing the most likely configuration for CRF
#' 
#' Approximate decoding with the block iterated conditional modes algorithm 
#' 
#' @param crf The CRF
#' @param blocks A list of vectors, each vector containing the nodes in a block
#' @param decode.method The decoding method to solve the clamped CRF
#' @param restart Non-negative integer to control how many restart iterations are repeated
#' @param start An initial configuration, a good start will significantly reduce the seraching time
#' @param ... The parameters for \code{decode.method}
#' @return This function will return the most likely configuration, which is a vector of length \code{crf$n.nodes}.
#' 
#' @examples
#' 
#' library(CRF)
#' data(Small)
#' d <- decode.block(Small$crf, list(c(1,3), c(2,4)))
#'
#' @export 
decode.block <- function(crf, blocks, decode.method = decode.tree, restart = 0, start = apply(crf$node.pot, 1, which.max), ...)
{
	y <- integer(crf$n.nodes)
	y[] <- start

	newcrf <- list()
	for (i in 1:length(blocks))
	{
		blocks[[i]] <- sort(blocks[[i]])
		clamped <- y
		clamped[blocks[[i]]] <- 0
		newcrf[[i]] <- clamp.crf(crf, clamped)
	}

	maxPot <- -1
	decode <- y
	restart <- max(0, restart)
	for (iter in 0:restart)
	{
		done = F
		while(!done)
		{
			done = T
			for (i in 1:length(blocks))
			{
				clamped <- y
				clamped[blocks[[i]]] <- 0
				newcrf[[i]] <- clamp.reset(newcrf[[i]], clamped)
				temp <- decode.method(newcrf[[i]], ...)
				if (any(temp != y[blocks[[i]]]))
				{
					y[blocks[[i]]] <- temp
					done = F
				}
			}
		}

		pot <- get.potential(crf, y)
		if (pot > maxPot)
		{
			maxPot <- pot
			decode <- y
		}

		if (iter < restart)
		{
			y <- ceiling(runif(crf$n.nodes) * crf$n.states)
		}
	}
	decode
}



#' Decoding method using integer linear programming
#' 
#' Computing the most likely configuration for CRF
#' 
#' Exact decoding with an integer linear programming formulation and approximate using LP relaxation
#' 
#' @param crf The CRF
#' @param lp.rounding Boolean variable to indicate whether LP rounding is need.
#' @return This function will return the most likely configuration, which is a vector of length \code{crf$n.nodes}.
#' 
#' @examples
#' 
#' library(CRF)
#' data(Small)
#' d <- decode.ilp(Small$crf)
#'
#' @export 
decode.ilp <- function(crf, lp.rounding = FALSE)
{
	vmap.nodes <- matrix(nrow=crf$n.nodes, ncol=2)
	vmap.edges <- matrix(nrow=crf$n.edges, ncol=2)
	n <- 0
	for (i in 1:crf$n.nodes)
	{
		vmap.nodes[i, 1] <- n + 1
		n <- n + crf$n.states[i]
		vmap.nodes[i, 2] <- n
	}
	vnum.nodes <- n
	for (i in 1:crf$n.edges)
	{
		vmap.edges[i, 1] <- n + 1
		n <- n + crf$n.states[crf$edges[i,1]] * crf$n.states[crf$edges[i,2]]
		vmap.edges[i, 2] <- n
	}
	vnum.edges <- n - vnum.nodes
	vnum.total <- n

	obj <- rep(0, vnum.total)
	for (i in 1:crf$n.nodes)
	{
		obj[vmap.nodes[i,1]:vmap.nodes[i,2]] <- -log(crf$node.pot[i,1:crf$n.states[i]])
	}
	for (i in 1:crf$n.edges)
	{
		obj[vmap.edges[i,1]:vmap.edges[i,2]] <- -log(crf$edge.pot[[i]])
	}
	obj[is.infinite(obj)] <- 1000

	cnum.nodes <- crf$n.nodes
	cnum.edges <- sum(crf$n.states[crf$edges])
	cnum.total <- cnum.nodes + cnum.edges
	mat <- matrix(0, nrow=cnum.total, ncol=vnum.total)
	for (i in 1:crf$n.nodes)
	{
		mat[i, vmap.nodes[i,1]:vmap.nodes[i,2]] <- 1
	}
	n <- cnum.nodes
	for (i in 1:crf$n.edges)
	{
		n1 <- crf$edges[i,1]
		n2 <- crf$edges[i,2]
		for (j in 1:crf$n.states[n1])
		{
			n <- n + 1
			mat[n, vmap.nodes[n1,1]-1+j] <- -1
			mat[n, seq.int(vmap.edges[i,1]-1+j, vmap.edges[i,2], by=crf$n.states[n1])] <- 1
		}
		for (j in 1:crf$n.states[n2])
		{
			n <- n + 1
			mat[n, vmap.nodes[n2,1]-1+j] <- -1
			mat[n, (vmap.edges[i,1]+(j-1)*crf$n.states[n1]):(vmap.edges[i,1]-1+j*crf$n.states[n1])] <- 1
		}
	}

	dir <- rep("==", cnum.total)
	rhs <- rep(0, cnum.total)
	rhs[1:cnum.nodes] <- 1

	if (lp.rounding)
	{
		types <- rep("C", vnum.total)
		bounds <- list(upper = list(ind = 1:vnum.total, val = rep(1, vnum.total)))
	}
	else
	{
		types <- rep("B", vnum.total)
		bounds <- NULL
	}

	result <- Rglpk::Rglpk_solve_LP(obj, mat, dir, rhs, types = types, bounds = bounds)

	if (result$status != 0)
	{
		warning("LP solution is not optimal.")
	}

	node.bel <- matrix(0, nrow=crf$n.nodes, ncol=crf$max.state)
	for (i in 1:crf$n.nodes)
	{
		node.bel[i, 1:crf$n.states[i]] <- result$solution[vmap.nodes[i,1]:vmap.nodes[i,2]]
	}
	apply(node.bel, 1, which.max)
}
