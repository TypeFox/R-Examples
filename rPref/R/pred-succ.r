
#' Predecessor and Successor Functions
#' 
#' Function for traversing the BTG (Better-Than-Graph or Hasse diagram) of a preference. 
#' 
#' @name pred_succ
#' @param df A data frame characterizing the set wherein predecessors/successors are searched.
#' @param v A numeric vector of indices in \code{df}. This represents the set of tuples for which predecessors/successors are searched.
#' @param p A preference. Worse tuples in the induced order are successors and better tuples are predecessors.
#' @param intersect Logical value. If it is \code{FALSE} (by default) the union of all predecessors/successors of \code{v} is returned.
#'        For \code{intersect = TRUE} the intersection of those values is returned.
#' 
#' @details
#' 
#' These functions return the predecessors and successors in the Better-Than-Graph of a preference.
#' Note that the successors/predecessors can can be plotted via \code{\link{get_btg}}. 
#' Before any of the successor/predecessor functions can be used the initialization has to be called as follows:
#' 
#' \code{init_pred_succ(p, df)}
#' 
#' There \code{p} is a preference object and \code{df} a data frame. 
#' This statement calculates the Better-Than-Relation on \code{df} w.r.t. \code{p}. 
#' Afterwards the subsequent predecessor and successor functions can be called. 
#' The value of \code{v} is a numeric vector within \code{1:nrow(df)} 
#' and characterizes a subset of tuples in \code{df}. 
#' The return value of these functions is again a numeric vector referring to the row numbers in \code{df} 
#' and it is always ordered ascending, independently of the order of the indices in \code{v}.
#' 
#' 
#' \describe{
#'   \item{\code{all_pred(p, v)}}{Returns all predecessors of \code{v}, i.e., indices of better tuples than \code{v}.}
#'   \item{\code{all_succ(p, v)}}{Returns all successors of \code{v}, i.e., indices of worse tuples than \code{v}.}
#'   \item{\code{hasse_pred(p, v)}}{Returns the direct predecessors of \code{v}, 
#'         i.e., indices of better tuples than \code{v} where the better-than-relation is contained in the transitive reduction.}
#'   \item{\code{hasse_succ(p, v)}}{Returns the direct successors of \code{v}, 
#'         i.e., indices of worse tuples than \code{v} where the better-than-relation is contained in the transitive reduction.}
#' }
#' 
#' If \code{v} has length 1, then the value of \code{intersect} does not matter, as there is nothing to intersect or join. 
#' For scalar values \code{x} and \code{y} the following identities hold, where \code{f} is one of the predecessor/successor functions:
#' 
#' \code{f(p, c(x, y), intersect = FALSE) == union(f(p, x), f(p, y))}
#' 
#' \code{f(p, c(x, y), intersect = TRUE) == intersect(f(p, x), f(p, y))}
#' 
#' 
#' @examples
#' 
#' # preference on mtcars for high mpg and low weight
#' p <- high(mpg) * low(wt)
#' init_pred_succ(mtcars, p)
#' 
#' # helper to show mpg/hp values
#' show_vals <- function(x) mtcars[x,c('mpg','wt')]
#' 
#' # pick some tuple "in the middle"
#' show_vals(10)
#' 
#' # show (direct) predecessors/successors of tuple 10
#' show_vals(hasse_pred(p, 10)) # Next better car
#' show_vals(hasse_succ(p, 10)) # Next worse car
#' show_vals(all_pred(p, 10))   # All better cars
#' show_vals(all_succ(p, 10))   # All worse cars
NULL


#' @rdname pred_succ
#' @export
init_pred_succ <- function(df, p) {
  p$init_pred_succ(df)
}

#' @rdname pred_succ
#' @export
hasse_pred <- function(p, v, intersect = FALSE) {
  p$h_predsucc(v, do_intersect = intersect, succ = FALSE)
}

#' @rdname pred_succ
#' @export
hasse_succ <- function(p, v, intersect = FALSE) {
  p$h_predsucc(v, do_intersect = intersect, succ = TRUE)
}

#' @rdname pred_succ
#' @export
all_pred <- function(p, v, intersect = FALSE) {
  p$all_predsucc(v, do_intersect = intersect, succ = FALSE)
}

#' @rdname pred_succ
#' @export
all_succ <- function(p, v, intersect = FALSE) {
  p$all_predsucc(v, do_intersect = intersect, succ = TRUE)
}