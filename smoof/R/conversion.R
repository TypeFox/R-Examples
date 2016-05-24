#' @title
#' Conversion between minimization and maximization problems.
#'
#' @description
#' We can minimize f by maximizing -f. The majority of predefined objective functions
#' in \pkg{smoof} should be minimized by default. However, there is a handful of
#' functions, e.g., Keane or Alpine02, which shall be maximized by default.
#' For benchmarking studies it might be beneficial to inverse the direction.
#' The functions \code{convertToMaximization} and \code{convertToMinimization}
#' do exactly that keeping the attributes.
#'
#' @note
#' Internally no wrapper is put around the original function. Instead the function
#' is copied and the body of the function is manipulated via the \code{body} function.
#' Both functions will quit with an error if multi-objective functions are passed.
#'
#' @param fn [\code{smoof_function}]\cr
#'  Smoof function.
#' @return [\code{smoof_function}]
#' @examples
#' # create a function which should be minimized by default
#' fn = makeSphereFunction(1L)
#' print(shouldBeMinimized(fn))
#' # Now invert the objective direction ...
#' fn2 = convertToMaximization(fn)
#' # and invert it again
#' fn3 = convertToMinimization(fn2)
#' # Now to convince ourselves we render some plots
#' opar = par(mfrow = c(1, 3))
#' plot(fn)
#' plot(fn2)
#' plot(fn3)
#' par(opar)
#' @name conversion
#' @rdname conversion
#' @export
convertToMaximization = function(fn) {
  convertProblemDirection(fn, minimize.after = FALSE)
}

#' @rdname conversion
#' @export
convertToMinimization = function(fn) {
  convertProblemDirection(fn, minimize.after = TRUE)
}

convertProblemDirection = function(fn, minimize.after = TRUE) {
  if (isMultiobjective(fn)) {
    stopf("Conversion to maximization only supported for single-objective problems
      at the moment, but your function '%s' has %i", getName(fn), getNumberOfObjectives(fn))
  }

  # If both are true, we want to convert min to min
  # If both are false, we want to convert max to max
  # Otherwise the conversion is ok
  if ((shouldBeMinimized(fn) && minimize.after) || (!shouldBeMinimized(fn) && !minimize.after)) {
    stopf("Function should already be %s.", (if (minimize.after) "minimized" else "maximized"))
  }

  # make copy of function
  fn2 = fn

  # multiply function body with -1
  body(fn2) = substitute(-1 * FUN, list(FUN = body(fn)))

  # copy attributes (get dropped on body() call)
  attributes(fn2) = attributes(fn)

  # flip maximization stuff
  fn2 = setAttribute(fn2, "minimize", !shouldBeMinimized(fn))
  return(fn2)
}
