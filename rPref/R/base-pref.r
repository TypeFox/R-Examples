#' Base Preferences
#' 
#' Base preferences are used to describe the different goals (dimensions, in case of a Skyline query)
#' of a preference query. 
#' 
#'
#' @name base_pref
#' @param expr A numerical/logical expression which is the term to evaluate for the current preference. 
#'       The objective is to search for minimal/maximal values of this expression (for \code{low}/\code{high}) or for 
#'       logical \code{TRUE} values (for \code{true}).
#' @param df (optional) A data frame, having the same structure (i.e., columns)
#'        like that data frame, where this preference is evaluated later on. 
#'        Causes a partial evaluation of the preference. Only the column names of \code{df} are relevant. 
#'        See below for details.
#' @param x An object to be tested if it is a base preference.
#' 
#' @details
#' 
#' Mathematically, all base preferences are strict weak orders (irreflexive, transitive and negative transitive).
#' 
#' The three fundamental base preferences are:
#' 
#' \describe{
#'   \item{\code{low(a), high(a)}}{Search for minimal/maximal values of \code{a}, 
#'         i.e., the induced order is the "smaller than" or "greater than" order on the values of \code{a}.
#'         The values of \code{a} must be numeric values.}
#'   \item{\code{true(a)}}{Searches for true values in logical expressions, i.e., \code{TRUE} is considered to be better than \code{FALSE}.
#'         The values of \code{a} must be logical values.
#'         For a tuplewise evaluation of a complex logical expression one has to use the \code{&} and \code{|} operators for logical AND/OR
#'         (and not the \code{&&} and \code{||} operators).}
#' }
#' 
#' 
#' The term \code{expr} may be just a single attribute or may contain an arbitrary expression,
#' depending on more than one attribute, e.g., \code{low(a+2*b+f(c))}.
#' There \code{a}, \code{b} and \code{c} are columns of the addressed data set and \code{f} has to be a previously defined function.
#' 
#' Functions contained in \code{expr} are evaluated over the entire data set, i.e., 
#' it is possible to use aggregate functions (\code{min}, \code{mean}, etc.). 
#' Note that all functions (and also variables which are not columns of the data set, where \code{expr} will be evaluated on)
#' must be defined in the same environment (e.g., environment of a function or global environment) as the base preference is defined.
#' 
#' The function \code{is.base_pref} returns \code{TRUE} if \code{x} is a preference object and \code{FALSE} otherwise.
#' 
#' 
#' @section Partial Evaluation of Preferences:
#' 
#' If the optional parameter \code{df} is given, 
#' then the expression is evaluated at the time of definition as far as possible.
#' All variables occurring as columns in \code{df} remain untouched. For example, consider
#' 
#' \code{f <- function(x) 2*x} \cr
#' \code{p <- true(cyl == f(1), mtcars)}
#' 
#' Then \code{p} is equivalent to the preference \code{true(cyl == 2)} as the variable \code{cyl} is a column in \code{mtcars}.
#' The rows of \code{df} are not relevant, e.g., using \code{mtcars[0,]} instead of \code{mtcars} makes no difference.
#' 
#' The preference selection, i.e., \code{psel(mtcars, p)} can be invoked without the partial evaluation.
#' But this results in an error, if the function \code{f} has meanwhile removed from the current environment.
#' Hence it is safer to do an early partial evaluation of all preferences, as far as they contain user defined functions.
#' 
#' The partial evaluation can be done manually by \code{\link{eval.pref}}.
#' 
#' 
#' @section Using Expressions in Preferences:
#' 
#' The \code{low_}, \code{high_} and \code{true_} preferences have the same functionality
#' as \code{low}, \code{high} and \code{true} 
#' but expect an expression \code{e} or symbol \code{e} as argument.
#' For example, \code{low(a)} is equivalent to \code{low_(expression(a))} or \code{low_(as.symbol("a"))}. 
#' 
#' 
#' This is very helpful for developing your own base preferences. Assume you want to define a base Preference \code{false}
#' as the dual of \code{true}. A definition like \code{false <- function(x) -true(x)} is the wrong approach, as 
#' \code{psel(data.frame(a = c(1,2)), false(a == 1))} will result in the error "object 'a' not found".
#' This is because \code{a} is considered as a variable and not as an (abstract) symbol to be evaluated later.
#' By defining
#' 
#' \code{false <- function(x, ...) -true_(substitute(x), ...)}
#' 
#' one gets a preference which behaves like a "built-in" preference.  
#' Additional optional parameters (like \code{df}) are bypassed.
#' The object \code{false(a == 1)} will output 
#' \code{[Preference] -true(a == 1)} on the console and 
#' \code{psel(data.frame(a = c(1,2)), false(a==1))} returns correctly the second tuple with \code{a==2}.
#' 
#' There is a special symbol \code{df__} which can be used in preference expression to access the given 
#' data set \code{df}, when \code{\link{psel}} is called on this data set. 
#' For example, on a data set where the first column has the name \code{A}
#' the preference \code{low(df__[[1]])} is equivalent to \code{low(A)}.
#' 
#' 
#' @seealso See \code{\link{complex_pref}} how to compose complex preferences to retrieve e.g., the Skyline.
#' See \code{\link{general_pref}} for functions applying to all kind of preferences.
#' See \code{\link{base_pref_macros}} for more base preferences.
#' 
#' @examples
#' # define a preference with a score value combining mpg and hp
#' p1 <- high(4 * mpg + hp)
#' # Perform the preference selection
#' psel(mtcars, p1)
#' 
#' # define a preference with a given function
#' f <- function(x, y) (abs(x - mean(x))/max(x) + abs(y - mean(y))/max(y))
#' p2 <- low(f(mpg, hp))
#' psel(mtcars, p2)
NULL


#' @rdname base_pref
#' @export
low <- function(expr, df = NULL) {
  expr <- as.expression(substitute(expr))
  return(eval.pref.internal(lowpref(expr, parent.frame()), df))
}

#' @rdname base_pref
#' @export
low_ <- function(expr, df = NULL) {
  return(eval.pref.internal(lowpref(as.expression(expr), parent.frame()), df))
}

#' @rdname base_pref
#' @export
high <- function(expr, df = NULL) {
  expr <- as.expression(substitute(expr))
  return(eval.pref.internal(highpref(expr, parent.frame()), df))
}

#' @rdname base_pref
#' @export
high_ <- function(expr, df = NULL) {
  return(eval.pref.internal(highpref(as.expression(expr), parent.frame()), df))
}

#' @rdname base_pref
#' @export
true <- function(expr, df = NULL) {
  expr <- as.expression(substitute(expr))
  return(eval.pref.internal(truepref(expr, parent.frame()), df))
}

#' @rdname base_pref
#' @export
true_ <- function(expr, df = NULL) {
  return(eval.pref.internal(truepref(as.expression(expr), parent.frame()), df))
}

#' @rdname base_pref
#' @export
is.base_pref <- function(x) {
  return(inherits(x, "basepref"))
}
