#' Make CRF
#' 
#' Generate CRF from the adjacent matrix
#' 
#' The function will generate an empty CRF from a given adjacent
#' matrix. If the length of \code{nstates} is less than \code{n.nodes}, it will
#' be used repeatly. All node and edge potentials are initilized as 1.
#' 
#' Since the CRF data are often very huge, CRF is implemented as an environment.
#' The assignment of environments will only copy the addresses instead of real data,
#' therefore the variables using normal assignment will refer to the exactly same CRF.
#' For complete duplication of the data, please use \code{\link{duplicate.crf}}. 
#' 
#' @param adj.matrix The adjacent matrix of CRF network.
#' @param n.states The state numbers of nodes.
#' @param n.nodes The number of nodes, which is only used to generate linear chain CRF when \code{adj.matrix} is NULL.
#' @return The function will return a new CRF, which is an environment with
#' components: 
#'   \item{n.nodes}{The number of nodes.} 
#'   \item{n.edges}{The number of edges.} 
#'   \item{n.states}{The number of states for each node. It is a vector of length \code{n.nodes}.} 
#'   \item{max.state}{The maximum number of states. It is equal to \code{max(n.states)}.} 
#'   \item{edges}{The node pair of each edge. It is a matrix with 2 columns and \code{n.edges} rows. Each row
#'     denotes one edge. The node with smaller id is put in the first column.}
#'   \item{n.adj}{The number of adjacent nodes for each node. It is a vector of length \code{n.nodes}.} 
#'   \item{adj.nodes}{The list of adjacent nodes for each
#'     node. It is a list of length \code{n.nodes} and the i-th element is a vector
#'     of length \code{n.adj[i]}.} 
#'   \item{adj.edges}{The list of adjacent edges for each node. It is similiar to \code{adj.nodes}
#'     while contains the edge ids instead of node ids.} 
#'   \item{node.pot}{The node potentials. It is a matrix with dimmension \code{(n.nodes, max.state)}.
#'     Each row \code{node.pot[i,]} denotes the node potentials of the i-th node.} 
#'   \item{edge.pot}{The edge potentials. It is a list of \code{n.edges} matrixes. Each matrix
#'     \code{edge.pot[[i]]}, with dimension \code{(n.states[edges[i,1]],
#'     n.states[edges[i,2]])}, denotes the edge potentials of the i-th edge.}
#'
#' @seealso \code{\link{duplicate.crf}}, \code{\link{clamp.crf}}, \code{\link{sub.crf}}
#'
#' @examples
#' 
#' library(CRF)
#' 
#' nNodes <- 4
#' nStates <- 2
#' 
#' adj <- matrix(0, nrow=nNodes, ncol=nNodes)
#' for (i in 1:(nNodes-1))
#' {
#' 	adj[i,i+1] <- 1
#' 	adj[i+1,i] <- 1
#' }
#' 
#' crf <- make.crf(adj, nStates)
#' 
#' crf$node.pot[1,] <- c(1, 3)
#' crf$node.pot[2,] <- c(9, 1)
#' crf$node.pot[3,] <- c(1, 3)
#' crf$node.pot[4,] <- c(9, 1)
#' 
#' for (i in 1:crf$n.edges)
#' {
#'    crf$edge.pot[[i]][1,] <- c(2, 1)
#'    crf$edge.pot[[i]][2,] <- c(1, 2)
#' }
#' 
#' @import Matrix
#' 
#' @export
make.crf <- function(adj.matrix = NULL, n.states = 2, n.nodes = 2)
{
	data <- new.env()
  
  if (is.null(adj.matrix))
    adj.matrix <- sparseMatrix(1:(n.nodes-1), 2:n.nodes, x = T, dims = c(n.nodes, n.nodes))
  
	if (length(dim(adj.matrix)) != 2 || dim(adj.matrix)[1] != dim(adj.matrix)[2])
		stop("Parameter 'adj.matrix' should be a square matrix")
	data$n.nodes <- dim(adj.matrix)[1]

	e <- which(adj.matrix != 0, arr.ind = TRUE)
	e <- matrix(c(e, e[,2], e[,1]), ncol=2)
	e <- unique(matrix(e[e[,1] < e[,2],], ncol=2))
	data$edges <- matrix(e[order(e[,1], e[,2]),], ncol=2)
	data$n.edges <- nrow(data$edges)

	.Call(Make_AdjInfo, data)

	data$n.states <- rep(n.states, length.out=data$n.nodes)
	data$max.state <- max(n.states)

	data$node.pot <- array(1, dim=c(data$n.nodes, data$max.state))
	data$edge.pot <- lapply(1:data$n.edges, function(i) array(1, dim=c(data$n.states[data$edges[i,1]], data$n.states[data$edges[i,2]])))

	class(data) <- "CRF"
	data
}



#' Duplicate CRF
#' 
#' Duplicate an existing CRF
#' 
#' This function will duplicate an existing CRF. Since CRF is implemented as an
#' environment, normal assignment will only copy the pointer instead of the 
#' real data. This function will generate a new CRF and really copy all data.
#' 
#' @param crf The existing CRF
#' @return The function will return a new CRF with copied data
#'
#' @seealso \code{\link{make.crf}}
#'
#' @export
duplicate.crf <- function(crf)
{
	data <- new.env()
	for (i in ls(envir=crf)) assign(i, get(i, envir=crf), envir=data)
	data
}



#' Calculate the potential of CRF
#' 
#' Calculate the potential of a CRF with given configuration
#' 
#' The function will calculate the potential of a CRF with given configuration,
#' i.e., the assigned states of nodes in the CRF.
#' 
#' @param crf The CRF
#' @param configuration The vector of states of nodes
#' @return The function will return the potential of CRF with given configuration
#'
#' @seealso \code{\link{get.logPotential}}
#'
#' @export
get.potential <- function(crf, configuration)
	.Call(Get_Potential, crf, configuration)



#' Calculate the log-potential of CRF
#' 
#' Calculate the logarithmic potential of a CRF with given configuration
#' 
#' The function will calculate the logarithmic potential of a CRF with given configuration,
#' i.e., the assigned states of nodes in the CRF.
#' 
#' @param crf The CRF
#' @param configuration The vector of states of nodes
#' @return The function will return the log-potential of CRF with given configuration
#'
#' @seealso \code{\link{get.potential}}
#'
#' @export
get.logPotential <- function(crf, configuration)
	.Call(Get_LogPotential, crf, configuration)
