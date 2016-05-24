## This file is part of the 'agop' library.
##
## Copyright 2013 Marek Gagolewski, Anna Cena
##
## 'agop' is free software: you can redistribute it and/or modify
## it under the terms of the GNU Lesser General Public License as published by
## the Free Software Foundation, either version 3 of the License, or
## (at your option) any later version.
##
## 'agop' is distributed in the hope that it will be useful,
## but WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
## GNU Lesser General Public License for more details.
##
## You should have received a copy of the GNU Lesser General Public License
## along with 'agop'. If not, see <http://www.gnu.org/licenses/>.


#' @title
#' Weak Dominance Relation (Preorder)
#' 
#' @description
#' \eqn{x <= y} iff \eqn{nx <= ny} AND
#' for all \eqn{i = 1,...,n} \eqn{x_{(n-i+1)} <= y_{(m-i+1)}},
#' where \eqn{nx = length(x)} and \eqn{ny = length(y)}.
#' 
#' @details
#' This function accepts only vectors with nonnegative elements.
#' 
#' @param x numeric vector
#' @param y numeric vector
#' @return single logical value; whether x <= y
#' 
#' @family binary_relations
#' @export
pord_weakdom <- function(x, y)
{
   .Call("pord_weakdom", x, y, PACKAGE="agop")
}


#' @title
#' Create Adjacency Matrix of Given Binary Relation
#' 
#' @description
#' Note that adjacency matrix  can also be conceived as a directed graph (DAG).
#' ret[i,j] iff i <= j
#' 
#' @param x list with elements to compare, preferrably named
#' @param pord function with 2 arguments, returning boolean value
#' @param ... additional arguments passed to \code{pord}
#' @return square 0-1 Matrix (of class \code{Matrix})
#' @family binary_relations
#' @export
rel_graph <- function(x, pord, ...)
{
   stopifnot(is.list(x))
   n <- length(x)
   ord <- matrix(0L, nrow=n, ncol=n)
   colnames(ord) <- names(x)
   rownames(ord) <- names(x)
   
   for (i in seq_along(x)) {
      for (j in seq_along(x)) {
         if (pord(x[[i]], x[[j]], ...))
            ord[i,j] <- 1L
      }
   }
   
   Matrix(ord)
}


#' @title
#' Check if Given Adjacency Matrix is Reflexive
#' 
#' @description
#' A binary relation R is reflexive, iff
#' for all \code{x} we have \eqn{xRx}.
#' The function just checks whether all elements
#' on the diagonal of \code{B} are non-zeros.
#' 
#' @param B object of class \code{igraph} or a square
#' 0-1 matrix of class \code{Matrix} or \code{matrix}
#' @return single logical value
#' @family binary_relations
#' @export
is_reflexive <- function(B)
{
   if (is(B, 'igraph')) B <- get.adjacency(B)
   stopifnot(is.matrix(B) || is(B, 'Matrix'), nrow(B) == ncol(B), nrow(B) > 0)
   all(Matrix::diag(B) != 0)
}


#' @title
#' Check if Given Adjacency Matrix is Total
#' 
#' @description
#' A binary relation R is total, iff
#' for all \code{x}, \code{y} we have \eqn{xRy} or \eqn{yRx}.
#' 
#' @param B object of class \code{igraph} or a square
#' 0-1 matrix of class \code{Matrix} or \code{matrix}
#' @return single logical value
#' @family binary_relations
#' @export
is_total <- function(B)
{
   if (is(B, 'igraph')) B <- get.adjacency(B)
   stopifnot(is.matrix(B) || is(B, 'Matrix'), nrow(B) == ncol(B), nrow(B) > 0)
   all(B + Matrix::t(B) != 0)
}



#' @title
#' Check if Given Adjacency Matrix is Transitive
#' 
#' @description
#' A binary relation R is transitive, iff
#' for all \code{x}, \code{y}, \code{z} we have \eqn{xRy} and \eqn{yRz}
#' => \eqn{xRz}
#' 
#' @param B object of class \code{igraph} or a square
#' 0-1 matrix of class \code{Matrix} or \code{matrix}
#' @return single logical value
#' @family binary_relations
#' @export
is_transitive <- function(B)
{
#    # version 0.1
#    if (is(B, 'igraph')) B <- get.adjacency(B)
#    stopifnot(is.matrix(B) || is(B, 'Matrix'), nrow(B) == ncol(B), nrow(B) > 0)
#    # slow as hell!
#    # @TODO - make faster
#    n <- nrow(B)
#    for (i in 1:n) {
#       for (j in 1:n) {
#          for (k in 1:n) {
#             if (as.logical(B[i,j]) &&  as.logical(B[j,k]) && !as.logical(B[i,k]))
#                return(FALSE)
#          }
#       }
#    }
#    TRUE
   
   # version 0.2
   if (!is(B, 'igraph')) B <- graph.adjacency(B)
   n <- vcount(B)
   for (i in 1:n) { # for each vertex
      # do breadth-first search from each vertex
      # transitivity holds iff each reachable vertex is within distance of 1
      disti <- graph.bfs(B, root=i, unreachable=FALSE, order=FALSE, dist=TRUE)$dist
      disti <- disti[!is.nan(disti)]
      if (length(disti) > 0 && any(disti > 1))
         return(FALSE)
   }
   TRUE # reached here -> no FALSE -> is transitive :-)
   
   
   # version 0.3 - todo
   # use shortest.paths()
}


#' @title
#' De-transitivitize Graph 
#' 
#' @description
#' Useful for draving Hasse diagrams.
#' 
#' 
#' @param B object of class \code{igraph} or a square
#' 0-1 matrix of class \code{Matrix} or \code{matrix}
#' @return object of class \code{Matrix}
#' @family binary_relations
#' @export
de_transitive <- function(B)
{
   if (is(B, 'igraph')) B <- get.adjacency(B)
   stopifnot(is.matrix(B) || is(B, 'Matrix'), nrow(B) == ncol(B), nrow(B) > 0)
   
   # slow as hell!
   n <- nrow(B)
   Matrix::diag(B) <- 0
   for (i in 1:n) {
      for (j in 1:n) {
         for (k in 1:n) {
            if (as.logical(B[i,j]) &&  as.logical(B[j,k]) && as.logical(B[i,k])
               && !as.logical(B[k,i]) && !as.logical(B[k,j]) && !as.logical(B[j,i]))
               B[i,k] <- 0
         }
      }
   }
   B
}



#' @title
#' Transitive Closure of Adjacency Matrix
#' 
#' @description
#' This may be slow for large graphs.
#' 
#' @param B object of class \code{igraph} or a square
#' 0-1 matrix of class \code{Matrix} or \code{matrix}
#' @return object of class \code{Matrix}
#' @export
#' @family binary_relations
closure_transitive <- function(B)
{
   if (is_transitive(B)) return(B) # do nothing
   if (is(B, 'igraph')) B <- get.adjacency(B)
   stopifnot(is.matrix(B) || is(B, 'Matrix'), nrow(B) == ncol(B), nrow(B) > 0)
   
   # slow as hell!
   n <- nrow(B)
   for (i in 1:n) {
      for (j in 1:n) {
         for (k in 1:n) {
            if (as.logical(B[i,j]) &&  as.logical(B[j,k]) && !as.logical(B[i,k]))
               B[i,k] <- 1
         }
      }
   }
   B
}



#' @title
#' Total Closure of Adjacency Matrix [Fair Totalization]
#' 
#' @description
#' Fair totalization: for each pair (x,y) s.t.
#' not xRy and not xRy let from now on xRy and yRx
#' 
#' @details
#' If you want a total preorder, call \code{\link{closure_transitive}}.
#' 
#' @param B object of class \code{igraph} or a square
#' 0-1 matrix of class \code{Matrix} or \code{matrix}
#' @return object of class \code{Matrix}
#' @export
#' @family binary_relations
closure_total_fair <- function(B)
{
   if (is(B, 'igraph')) B <- get.adjacency(B)
   stopifnot(is.matrix(B) || is(B, 'Matrix'), nrow(B) == ncol(B), nrow(B) > 0)
   
   B2 <- B + Matrix::t(B)
   # 0s lie symmetricly over diagonal and indicate incomparable pairs
   wh <- Matrix::which(B2 == 0, arr.ind=TRUE)
   B[wh] <- 1
   B
}



#' @title
#' Get Incomparable Pairs in an Adjacency Matrix
#' 
#' @description
#' A pair (x,y) is incomparable iff
#' not xRy and not xRy
#'  
#' @details
#' See also \code{\link{get_independent_sets}} of how to generate
#' all maximal independent sets.
#' 
#' @param B object of class \code{igraph} or a square
#' 0-1 matrix of class \code{Matrix} or \code{matrix}
#' @return integer matrix with two columns (indices of incomparable elements,
#' not that these are pairs, and not sets: you'll get (i,j) and (j,i))
#' @export
#' @family binary_relations
get_incomparable_pairs <- function(B)
{
   if (is(B, 'igraph')) B <- get.adjacency(B)
   stopifnot(is.matrix(B) || is(B, 'Matrix'), nrow(B) == ncol(B), nrow(B) > 0)
   
   B2 <- B + Matrix::t(B)
   # 0s lie symmetricly over diagonal and indicate incomparable pairs
   Matrix::which(B2 == 0, arr.ind=TRUE)
}


#' @title
#' Get All Maximal Independent Sets
#' 
#' @description
#' The function generates vectors of indices \eqn{S_j=\{i_1,...,i_{k_j}\}}
#' such that all pairs from \eqn{S_j} are incomparable
#' (A pair (i,i') is incomparable iff
#' not \eqn{i R i'} and not \eqn{i' R i},
#' see also \code{\link{get_incomparable_pairs}}.
#' 
#' 
#' @details
#' Note that we assume that \eqn{B} is transitive.
#' Loops are not taken into account at all.
#' 
#' @param B object of class \code{igraph} or a square
#' 0-1 matrix of class \code{Matrix} or \code{matrix}
#' @return list of integer vectors; each list element defines
#' an independent set of vertices numbers
#' 
#' @export
#' @family binary_relations
get_independent_sets <- function(B)
{
   if (!is(B, 'igraph')) B <- graph.adjacency(B)
   
   # create graph with all not-directly-connected nodes
   C <- (shortest.paths(B, mode="all") > 1) # symmetric
   diag(C) <- 0
   
   # find all maximal cliques (== independent sets)
   out <- maximal.cliques(graph.adjacency(C, mode="undirected"))
   out[sapply(out, length) > 1] # remove singletons
}



#' @title
#' Get All Equivalence Classes of a Total Binary Relation
#' 
#' 
#' @description
#' Note that we assume that \eqn{B} is total, reflexive and transitive.
#' 
#' @param B object of class \code{igraph} or a square
#' 0-1 matrix of class \code{Matrix} or \code{matrix}
#' @return list of integer vectors; each list element defines
#' an equivalence class by listing vertices' numbers;
#' each vector is ordered by the outdegrees of their nodes (they are
#' the same in each class)
#' 
#' @export
#' @family binary_relations
get_equivalence_classes <- function(B)
{
   if (!is(B, 'igraph')) { BAD <- as.matrix(B); B <- graph.adjacency(B) }
   else { BAD <- as.matrix(get.adjacency(B)) }
   stopifnot(is_total(BAD), is_transitive(BAD))
   
   # now B - igraph
   # now BAD - matrix
   
   n <- nrow(BAD)
   
   # create graph with all `symmetric` pairs
   C <- BAD & t(BAD) # logical AND
   out <- maximal.cliques(graph.adjacency(C, mode="undirected"))
   
   # now sort according to the out-degree (decreasingly)
   # every vertex in an equivalence class has the same degree
   # if x<y, then outdegree(x)>outdegree(y)
   extractFirst <- sapply(out, function(e) e[1])
   ord <- order(rowSums(BAD)[extractFirst], decreasing=TRUE)
   out[ord]
}


