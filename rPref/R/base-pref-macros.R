#' Useful Base Preference Macros
#' 
#' In addition to the fundamental base preferences, 
#' rPref offers some macros to define preferences where a given interval or point is preferred. 
#'
#' @name base_pref_macros
#' @param expr A numerical expression (for \code{around} and \code{between})
#'        or an arbitrary expression (for \code{pos} and \code{layered}).
#'        The objective are those tuples where \code{expr} evaluates to a value within the preferred interval, layer, etc. 
#'        Regarding attributes, functions and variables, the same requirements as for \code{\link{base_pref}} apply.
#' @param center Preferred numerical value for \code{around}.
#' @param left Lower limit (numerical) of the preferred interval for \code{between}.
#' @param right Upper limit (numerical) of the preferred interval for \code{between}.
#' @param pos_value A vector containing the preferred values for a \code{pos} preference.
#'         It has to be of the same type (numeric, logical, character, ...) as \code{expr}.
#' @param ... Layers (sets) for a \code{layered} preference. Each parameter corresponds to a layer 
#'            and the first one characterizes the most preferred values.
#' @param df (optional) Data frame for partial evaluation. See \code{\link{base_pref}} for details.
#' 
#' @section Definition of the Preference Macros:
#' 
#' \describe{
#'   \item{\code{between(expr, left, right)}}{Those tuples are preferred where \code{expr} evaluates
#'     to a value between \code{left} and \code{right}.
#'     For values not in this interval, the values nearest to the interval are preferred.}
#'   \item{\code{around(expr, center)}}{Same as \code{between(expr, center, center)}.}
#'   \item{\code{pos(expr, pos_value)}}{Those tuples are preferred, where \code{expr} evaluates 
#'     to a value which is contained in \code{pos_value}.}
#'   \item{\code{layered(expr, layer1, layer2, ..., layerN)}}{For the most preferred tuples \code{expr}
#'     must evaluate to a value in \code{layer1}. 
#'     The second-best tuples are those where \code{expr} evaluates to a value in \code{layer2} and so forth. 
#'      Values occurring in none of the layers are considered worse than those in \code{layerN}.
#'      Technically, this is realized by a prioritization chain (lexicographical order)
#'      of \code{\link{true}} preferences.}
#' }
#' 
#' Note that only the argument \code{expr} may contain columns from the data frame, 
#' all other variables must evaluate to explicit values. 
#' For example \code{around(mpg, mean(mpg))} is not allowed. In this case, one can use 
#' \code{around(mpg, mean(mtcars$mpg))} instead. Or alternatively, without using the base preference macros, 
#' \code{low(abs(mpg - mean(mpg)))} does the same. There, the actual mean value of \code{mpg} is calculated 
#' just when the preference selection via \code{psel} is called.
#'
#' @examples 
#' # search for cars where mpg is near to 25
#' psel(mtcars, around(mpg, 25))
#' 
#' # cyl = 2 and cyl = 4 are equally good, cyl = 6 is worse
#' psel(mtcars, layered(cyl, c(2, 4), 6))
NULL


#' @rdname base_pref_macros
#' @export
around <- function(expr, center, df = NULL) {
  expr <- as.expression(call('abs', call('-', substitute(expr), center)))
  return(eval.pref.internal(lowpref(expr, parent.frame()), df))
}


#' @rdname base_pref_macros
#' @export
between <- function(expr, left, right, df = NULL) {
  expr <- substitute(expr)
  between_expr <- as.expression(call('pmax', call('-', left, expr), 0, call('-', expr, right)))
  return(eval.pref.internal(lowpref(between_expr, parent.frame()), df))
}

#' @rdname base_pref_macros
#' @export
pos <- function(expr, pos_value, df = NULL) {
  expr <- as.expression(call('%in%', substitute(expr), pos_value))
  return(eval.pref.internal(truepref(expr, parent.frame()), df))
}

#' @rdname base_pref_macros
#' @export
layered <- function(expr, ..., df = NULL) {
  sl <- list(...)
  if (length(sl) < 2) stop("Empty layered preference is not allowed!")
  expr <- substitute(expr)
  layers <- lapply(sl, function(x) { 
    tmp_expr <- as.expression(call('%in%', expr, x))
    return(eval.pref.internal(truepref(tmp_expr, parent.frame()), df))
  })
  return(Reduce('&', layers))
}
