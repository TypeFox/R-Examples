## list_utils.R
##   - Utility functions for lists and vectors
##
## RGP - a GP system for R
## 2010 Oliver Flasch (oliver.flasch@fh-koeln.de)
## with contributions of Thomas Bartz-Beielstein, Olaf Mersmann and Joerg Stork
## released under the GPL v2
##

##' Functions for Lisp-like list processing
##'
##' Simple wrapper functions that allow Lisp-like list processing in
##' R: \code{first} to \code{fifth} return the first to fifth element
##' of the list \code{x}. \code{rest} returns all but the first
##' element of the list \code{x}. \code{is.empty} returns \code{TRUE}
##' iff the list \code{x} is of length 0. \code{is.atom} returns
##' \code{TRUE} iff the list \code{x} is of length 1.
##' \code{is.composite} returns \code{TRUE} iff the list \code{x} is
##' of length > 1. \code{contains} return \code{TRUE} iff the list
##' \code{x} contains an element identical to \code{elt}.
##'
##' @param x A list or vector.
##' @param elt An element of a list or vector.
##'
##' @rdname lispLists
first <- function(x) x[[1]]

##' @rdname lispLists
rest <- function(x) x[-1]

##' @rdname lispLists
second <- function(x) x[[2]]

##' @rdname lispLists
third <- function(x) x[[3]]

##' @rdname lispLists
fourth <- function(x) x[[4]]

##' @rdname lispLists
fifth <- function(x) x[[5]]

##' @rdname lispLists
is.empty <- function(x) length(x) == 0

##' @rdname lispLists
is.atom <- function(x) length(x) == 1

##' @rdname lispLists
is.composite <- function(x) length(x) > 1

##' @rdname lispLists
contains <- function(x, elt) {
  for (candidate in x) if (identical(candidate, elt)) return(TRUE)
  FALSE
}

##' Sort a vector or list by the result of applying a function
##'
##' Sorts a vector or a list by the numerical result of applying the function \code{byFunc}.
##'
##' @param xs A vector or list.
##' @param byFunc A function from elements of \code{xs} to \code{numeric}.
##' @return The result of sorting \code{xs} by \code{byfunc}.
sortBy <- function(xs, byFunc) {
  if (missing(byFunc))
    o <- order(xs)
  else
    o <- order(sapply(xs, byFunc))
  xs[o]
}

##' Sort a vector or list via a given ranking
##'
##' Reorders a vector or list according to a given ranking \code{ranking}.
##'
##' @param xs The vector or list to reorder.
##' @param ranking The ranking to sort \code{xs} by, defaults to \code{rank(xs)}.
##' @return The result of reordering \code{xs} by \code{ranking}.
sortByRanking <- function(xs, ranking = rank(xs)) {
  sorted <- xs
  sorted[ranking] <- xs
  sorted
}

##' Sorting algorithms for vectors and lists
##'
##' These algorithms sort a list or vector by a given order relation (which
##' defaults to \code{<=}).
##' \code{insertionSort} is a stable O(n^2) sorting algorithm that is quite efficient
##' for very small sets (less than around 20 elements). Use an O(n*log(n)) algorithm
##' for larger sets.
##'
##' @param xs The vector or list to sort.
##' @param orderRelation The orderRelation to sort \code{xs} by (defaults to \code{`<=`}).
##'   This relation by should reflexive, antisymetric, and transitive.
##' @return The vector or list \code{xs} sorted by the order relation \code{orderRelation}.
##'
##' @rdname sortingAlgorithms
##' @export
insertionSort <- function(xs, orderRelation = NULL) {
  orderRelation <- if (is.null(orderRelation)) `<=` else orderRelation
  sorted <- xs
  l <- length(sorted)
  if (l <= 1) {
    sorted
  } else {
    for (i in 2:l) {
      xi <- xs[[i]]
      j <- i
      while (j > 1 &&  sorted[[j - 1]] > xi) {
        sorted[[j]] <- sorted[[j - 1]]
        j <- j - 1
      }
      sorted[[j]] <- xi
    }
    sorted
  }
}

##' Calculate the inverse of a permutation
##'
##' Returns the inverse of a permutation \code{x} given as an integer vector.
##' This function is useful to turn a ranking into an ordering and back, for example.
##'
##' @param x The permutation to return the inverse for.
##' @return The inverse of the permutation \code{x}.
##' @seealso \code{\link{rank}}, \code{\link{order}}
inversePermutation <- function(x) {
  l <- length(x); o <- numeric(l)
  o[x] <- 1:l
  o
}

##' Choose a random element from a list or vector
##'	
##' Returns a unformly random chosen element of the vector or list \code{x}.
##'
##' @param x The vector or list to chose an element from.
##' @param prob A vector of probability weights for obtaining the elements of the
##'   vector or list being sampled.
##' @return A uniformly random element of \code{x}.
randelt <- function(x, prob = NULL) {
  l <- length(x)
  if (l == 0)
    NULL
  else
    sample(x, 1, prob = prob)[[1]] # sample(x, 1) always yields a list of exactly one element
}

##' Splitting and grouping of lists
##'
##' Functions for splitting and grouping lists into sublists.
##' \code{splitList} splits a list \code{l} into \code{max(groupAssignment)} groups.
##' The integer indices of \code{groupAssignment} determine in which group each
##' element of \code{l} goes.
##' \code{groupListConsecutive} splits \code{l} into \code{numberOfGroups} consecutive
##' sublists (or groups).
##' \code{groupListDistributed} distributes \code{l} into \code{numberOfGroups}
##' sublists (or groups).
##' \code{flatten} flattens a list \code{l} of lists into a flat list by concatenation. If
##' \code{recursive} is \code{TRUE} (defaults to \code{FALSE}), flatten will be recursively
##' called on each argument first.
##' \code{intersperse} joins two lists \code{xs} and \code{ys} into a list of pairs containig
##' every possible pair, i.e. \code{intersperse(xs, ys)} equals the product list of \code{xs}
##' and \code{ys}. The \code{pairConstructor} parameter can be used to change the type of pairs
##' returned.
##'
##' @param l A list.
##' @param xs A list.
##' @param ys A list.
##' @param pairConstructor The function to use for constructing pairs, defaults to \code{list}.
##' @param groupAssignment A vector of group assignment indices.
##' @param numberOfGroups The number of groups to create, must be <= length(l)
##' @param recursive Whether to operate recursively on sublists or vectors.
##' @return A list of lists, where each member represents a group.
##'
##' @rdname listSplittingAndGrouping
splitList <- function(l, groupAssignment) {
  groups <- Map(function(i) list(), 1:max(groupAssignment))
  for (i in 1:length(groupAssignment))
    groups[[groupAssignment[i]]] <- c(groups[[groupAssignment[i]]], l[[i]])
  groups         
}

##' @rdname listSplittingAndGrouping
groupListConsecutive <- function(l, numberOfGroups) {
  ll <- length(l)
  if (numberOfGroups > ll) stop("groupListConsecutive: numberOfGroups must be <= length(l)")
  gl <- floor(ll / numberOfGroups)
  glr <- ll %% numberOfGroups
  gassign <- c(Reduce(c, lapply(1:numberOfGroups, function(gi) rep(gi, gl))), rep(numberOfGroups, glr))  
  splitList(l, gassign)       
}

##' @rdname listSplittingAndGrouping
groupListDistributed <- function(l, numberOfGroups) {
  ll <- length(l)
  if (numberOfGroups > ll) stop("groupListDistributed: numberOfGroups must be <= length(l)")
  gl <- floor(ll / numberOfGroups)
  gassign <- 1:ll
  for (i in 1:ll) gassign[i] <- (i - 1) %% numberOfGroups + 1
  splitList(l, gassign)
}

##' @rdname listSplittingAndGrouping
flatten <- function(l, recursive = FALSE)
  if (recursive && !is.atom(l)) {
    Reduce(function(a, b) c(flatten(a, recursive = TRUE), flatten(b, recursive = TRUE)), l, init = list()) 
  } else {
    Reduce(c, l, init = list())
  }

##' @rdname listSplittingAndGrouping
intersperse <- function(xs, ys, pairConstructor = list)
  Reduce(c, Map(function(x1) Map(function(x2) pairConstructor(x1, x2), ys), xs))
