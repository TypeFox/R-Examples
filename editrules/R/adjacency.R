
#' Derive adjecency matrix from collection of edits
#'
#' A set of edits can be represented as a graph where every vertex is
#' an edit. Two vertices are connected if they have at least one variable
#' in \code{vars} in common.
#' 
#' \code{adjacency} returns the adjacency matrix. The elements of the matrix
#' count the number of variables shared by the edits indicated in the row- and 
#' column names. The adjacency matrix can be converted to an igraph object with 
#' \code{graph.adjacency}from the \code{igraph} package.
#'
#' \code{as.igraph} converts a set of edits to an \code{igraph} object directly.
#'
#'
#' @param E \code{\link{editmatrix}}, \code{\link{editarray}} or \code{\link{editset}}
#' @param nodetype adjacency between rules, vars or both?
#' @param rules selection of edits
#' @param vars selection of variables
#' @param ... arguments to be passed to or from other methods
#'
#' @return the adjacency matrix of edits in \code{E} with resect to 
#'  the variables in \code{vars}
#'
#' @example ../examples/graph.R
#'
#' @seealso \code{\link{plot.editmatrix}}, \code{\link{plot.editarray}}, \code{\link{plot.editset}}
#' @export  
adjacency <- function(E, nodetype=c("all", "rules","vars"), rules=rownames(E), vars=getVars(E),...){
    stopifnot( all(vars %in% getVars(E)) )
    UseMethod('adjacency')
}

#' @rdname adjacency
#' @method adjacency editmatrix
#' @export
adjacency.editmatrix <- function(E, nodetype=c("all", "rules","vars"), rules=rownames(E), vars=getVars(E),...){
    A <- abs(sign(getA(E)))
    nodetype <- match.arg(nodetype)
    adjec(A,nodetype=nodetype, rules=rules, vars=vars)
}


#' @rdname adjacency
#' @method adjacency editarray
#' @export
adjacency.editarray <- function(E, nodetype=c("all", "rules","vars"), rules=rownames(E), vars=getVars(E),...){
    A <- contains(E)
    nodetype <- match.arg(nodetype)
    adjec(A,nodetype=nodetype, rules=rules, vars=vars)
}


#' @rdname adjacency
#' @method adjacency editset
#' @export
adjacency.editset <- function(E, nodetype=c("all", "rules","vars"), rules=c(rownames(E$num),rownames(E$mixcat)), vars=getVars(E),...){
    A <- contains(E)
    nodetype <- match.arg(nodetype)
    adjec(A,nodetype=nodetype, rules=rules, vars=vars)
}





# derive adjacency from 1/0 or boolean matrix.
# Internal loops only, but nrow(A)^2 memory complexity. 
# future optimization options: sparse matrices, lower/upper triangle only.
adjec <- function(A, nodetype="all", rules=rownames(A), vars=colnames(A)){
  A <- A[rules, vars, drop=FALSE]
  m <- NULL
  vars <- NULL
  if (nodetype=="all"){
    N <- nrow(A) + ncol(A)
    nms <- c(rownames(A), colnames(A))
    vars <- rep(c(FALSE, TRUE), times=c(nrow(A), ncol(A)))
    m <- matrix(0, nrow=N, ncol=N, dimnames=list(nms, nms))
    m[!vars, vars] <- A
    m[vars, !vars] <- t(A)
  } else{
    vars <- rep(FALSE, nrow(A))
    if (nodetype=="vars"){
      vars <- rep(TRUE, ncol(A))
      A <- t(A)
    }
    I <- rep(1:nrow(A), times=nrow(A))
    J <- rep(1:nrow(A), each=nrow(A)) 
    m <- matrix(
            rowSums(A[I,,drop=FALSE] & A[J,,drop=FALSE]),
            nrow=nrow(A),
            dimnames=list(rownames(A), rownames=rownames(A))
    )
    diag(m) <- 0
  }
  attr(m,"vars") <- vars
  m
}
