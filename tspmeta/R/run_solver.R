#' Runs a solver on a TSP instance.
#'
#' Currently the following solvers are supported:
#'   nearest_insertion: See \code{\link[TSP]{solve_TSP}}.
#'   farthest_insertion : See \code{\link[TSP]{solve_TSP}}.
#'   cheapest_insertion : See \code{\link[TSP]{solve_TSP}}.
#'   arbitrary_insertion: See \code{\link[TSP]{solve_TSP}}.
#'   nn: See \code{\link[TSP]{solve_TSP}}.
#'   repetitive_nn: See \code{\link[TSP]{solve_TSP}}.
#'   concorde: See \code{\link[TSP]{solve_TSP}}.
#'
#' @param x [\code{\link{tsp_instance}}]\cr
#'   TSP instance.
#' @param method [\code{character(1)}]\cr
#'   Solver to use on TSP instance. To use concorde and/or linkern it is necessary to specify the path
#'   to the concorde/linkern executable with \code{\link[TSP]{concorde_path}}.
#' @param ... [any]\cr
#'   Control parameters for solver.
#' @return [\code{\link[TSP]{TOUR}}]
#'   TOUR object from package TSP, containing order of cities, tour length and
#'   method name that generated this solution.
#'
#' @examples
#' x = random_instance(10)
#' tours = sapply(c("nn", "cheapest_insertion", "arbitrary_insertion"), function(solver) {
#'   list(solver = run_solver(x, method = solver))
#' })
#' \dontrun{
#'   concorde_path(path = "/absolute/path/to/concorde/executable")
#'   concorde_tour = run_solver(x, method = "concorde")
#'   concorde_tour = run_solver(x, method = "linkern")
#' }
#'
#' @export
run_solver = function(x, method, ...) {
  assertClass(x, "tsp_instance")
  assertChoice(method, choices = get_solvers())
  if (method %in% c("nearest_insertion",
                    "farthest_insertion",
                    "cheapest_insertion",
                    "arbitrary_insertion",
                    "nn",
                    "repetitive_nn",
                    "concorde",
                    "linkern")) {
    solve_TSP(as_TSP(x), method = method, control = list(...))
  } else if(method == "2-opt") {
    fast_two_opt(x)
  }
}

#' Returns integrated solver names.
#'
#' @return [\code{character}].
#' @export
get_solvers = function() {
  c(
    "nearest_insertion",
    "farthest_insertion",
    "cheapest_insertion",
    "arbitrary_insertion",
    "nn",
    "repetitive_nn",
    "2-opt",
    "concorde",
    "linkern"
  )
}

#' Runs 2-Opt local search on TSP instance.
#'
#' @param x [\code{\link{tsp_instance}}]\cr
#'   TSP instance.
#' @param initial_tour [\code{numeric}]\cr
#'   Initial tour.
#' @return [\code{\link[TSP]{TOUR}}]
#'   TOUR object from package TSP, containing order of cities, tour length and
#'   method name that generated this solution.
#' @useDynLib tspmeta do_fast_two_opt
fast_two_opt = function(x, initial_tour) {
  dists = as.matrix(x$dists)

  if (missing(initial_tour))
    initial_tour = sample(1:number_of_cities(x))

  res = .Call(do_fast_two_opt, dists, initial_tour)
  TOUR(x = res[[1]], method = "2-opt", tsp = as_TSP(x))
}
