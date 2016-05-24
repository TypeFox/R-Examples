# Generates crossproduct of exhaustive options joined with rows of a given design data.frame.
# State-based object with iterator functions is returned.
# If both arguments are empty, the a design is created with corresponds to one 'point' that
# defines a function call with no arguments.
# @param ex [\code{list}]\cr
#   Named list of parameters settings which should be exhaustively tried.
#   All elements of the list must be primitive vectors like numeric, integer, factor, etc.
# @param .design [\code{data.frame}]\cr
#   The design. Rows define one 'point'.
# @return List of funs. nextElem returns a named list (ordered by list element names).
designIterator = function(ex, .design = data.frame()) {
  nextState = function(state, pos = 1L) {
    if (state[pos] < state.last[pos])
      replace(state, pos, state[pos] + 1L)
    else
      nextState(replace(state, pos, 1L), pos + 1L)
  }

  nextElem = function() {
    state <<- nextState(state)
    counter <<- counter + 1L

    x = c(as.list(.design[state[! is.ex.state], , drop = FALSE]),
          mapply(function(n, s) ex[[n]][s], n = names.ex.state, s = state[is.ex.state], SIMPLIFY = FALSE))
    x[order(names2(x))]
  }

  hasNext = function() {
    counter < counter.max
  }

  reset = function() {
    state <<- state.init
    counter <<- 0L
    invisible(TRUE)
  }


  state.last = sort(setNames(c(vapply(ex, length, 1L), max(nrow(.design), 1L)), c(names(ex), ".design.row")), decreasing = TRUE)
  state.init = setNames(c(0L, rep.int(1L, length(state.last) - 1L)), names(state.last))
  counter.max = prod(state.last)
  if (counter.max > .Machine$integer.max)
    stop("The generated design is too big. Designs with up to ",
         .Machine$integer.max, " rows are supported!")
  counter.max = as.integer(counter.max)
  is.ex.state = (names(state.init) != ".design.row")
  names.ex.state = names(state.init)[is.ex.state]

  state = state.init
  counter = 0L

  list(nextElem = nextElem,
       hasNext = hasNext,
       reset = reset,
       n.states = counter.max,
       storage = c(vapply(.design, storage.mode, character(1L)),
                   vapply(ex, storage.mode, character(1L))))
}

#' @title Create parameter designs for problems and algorithms.
#'
#' @description
#' Create a parameter design for either a problem or an algorithm that you
#' can use in \code{\link{addExperiments}}.
#' All parameters in \code{design} and \code{exhaustive} be \dQuote{primitive}
#' in the sense that either \code{is.atomic} is \code{TRUE} or \code{is.factor} is \code{TRUE}.
#'
#' Be aware of R's default behaviour of converting strings into factors if you use the \code{design}
#' parameter. See option \code{stringsAsFactors} in \code{\link{data.frame}} to turn this off.
#' @param id [\code{character(1)}]\cr
#'
#'   Id of algorithm or problem.
#' @param design [\code{data.frame}]\cr
#'   The design. Must have named columns corresponding to parameters.
#'   Default is an empty \code{data.frame()}.
#' @param exhaustive [\code{list}]\cr
#'   Named list of parameters settings which should be exhaustively tried.
#'   Names must correspond to parameters.
#'   Default is empty list.
#' @return [\code{\link{Design}}].
#' @export
#' @aliases Design
#' @examples \dontrun{
#' # simple design for algorithm "a1" with no parameters:
#' design = makeDesign("a1")
#'
#' # design for problem "p1" using predefined parameter combinations
#' design = makeDesign("p1", design = data.frame(alpha = 0:1, beta = c(0.1, 0.2)))
#'
#' # creating a list of designs for several algorithms at once, all using the same
#' # exhaustive grid of parameters
#' designs = lapply(c("a1", "a2", "a3"), makeDesign,
#'                   exhaustive = list(alpha = 0:1, gamma = 1:10/10))
#' }
makeDesign = function(id, design = data.frame(), exhaustive = list()) {
  # ... if we had the registry here, we could do some sanity checks, e.g.
  # test if the storage mode of parameters matches the storage mode of those
  # in the database
  # if we push out a not backward compatible version, do this here.
  assertString(id)
  assertDataFrame(design, types = "atomic")
  assertList(exhaustive, types = "atomic", names = "named")
  if (any(viapply(exhaustive, length) == 0L))
    stop("All elements of exhaustive must have at least have length 1!")
  if (anyDuplicated(c(names(design), names(exhaustive))) > 0L)
    stop("Duplicated design parameters found!")
  setClasses(list(id = id, designIter = designIterator(exhaustive, .design = design)), "Design")
}

#' @export
print.Design = function(x, ...) {
  n = x$designIter$n.states
  storage = x$designIter$storage
  catf("Design for %s with %i row%s", x$id, n, ifelse(n == 1L, "", "s"))
  cat(collapse(sprintf("  %-10s: %s", names(storage), storage), "\n"), "\n")
}
