
#' Complex Preferences
#' 
#' Complex preferences are used to compose different preference orders. 
#' For example the Pareto composition (via operator \code{*}) is the usual operator
#' to compose the preference for a Skyline query. The Skyline is also known as Pareto frontier.
#' All complex preferences are mathematically strict partial orders (irreflexive and transitive).
#' 
#' @name complex_pref
#' @param p,p1,p2,x Preferences (they can be either base preferences, see \code{\link{base_pref}}, or complex preferences),
#'                  or, for \code{is.complex_pref}, an object to be tested if it is a complex preference.
#'                
#' 
#' @section Skylines:
#' 
#' The most important preference composition operator is the Pareto operator (\code{p1 * p2}) to formulate a Skyline query. 
#' A tuple t1 is better than t2 w.r.t. \code{p1 * p2} if it is strictly better w.r.t. one of the preferences p1, p2 and is better or equal w.r.t. the other preference.
#' 
#' The syntactical correspondence to other query languages supporting Skylines/preferences to rPref
#' is given as follows:
#' 
#' \itemize{
#' \item A query in the syntax from Borzsonyi et. al (2001) like
#' 
#' "\code{... SKYLINE OF a MAX, b MIN, c MAX}" 
#' 
#' corresponds in rPref to the preference
#' 
#' \code{high(a) * low(b) * high(c)}.
#' 
#' 
#' \item A query in the syntax from Kiessling (2002) like 
#' 
#' "\code{... PREFERRING a LOWEST AND (b HIGHEST PRIOR TO c LOWEST)}" 
#' 
#' corresponds in rPref to 
#' 
#' \code{low(a) * (high(b) & low(c))}.
#' 
#' 
#' \item A query in the syntax of the "Skyline" feature of the commercial database "EXASOL EXASolution 5" like 
#' 
#' "\code{... PREFERRING LOW a PLUS (b = 1 PRIOR TO LOW c))}" 
#' 
#' corresponds in rPref to
#' 
#' \code{low(a) * (true(b == 1) & low(c))}.
#' 
#' }
#' 
#' Note that preferences in rPref can be translated to some of this query dialects by \code{\link{show.query}}. 
#' 
#' @section Definition of Additional Preference Operators:
#' 
#' Additionally, rPref supports the following preference composition operators:
#' 
#' \describe{
#'  \item{\code{p1 & p2}}{Prioritization (lexicographical order): A tuple t1 is better than t2 w.r.t. \code{p1 & p2} if it is 
#'   strictly better w.r.t. \code{p1} or is equal w.r.t. \code{p1} and is better w.r.t. \code{p2}.}
#'  \item{\code{p1 | p2}}{Intersection preference: A tuple t1 is better than t2 w.r.t. \code{p1 | p2} if it is strictly better w.r.t. both preferences. 
#'        This is a stricter variant of the Pareto operator. The evaluation of \code{psel(df, p1 | p2)} is always a subset of \code{psel(df, p1 * p2)}.}
#'  \item{\code{p1 + p2}}{Union preference: A tuple t1 is better than t2 w.r.t. \code{p1 + p2} if it is strictly better w.r.t. to one of the preferences. 
#'  Note that this can violate the strict partial order property, if the domains (the tuples on which \code{p1} and \code{p2} define better-than-relationships) 
#'  of the preferences are not disjoint.}
#'  \item{\code{reverse(p1)} or \code{-p1}}{Reverse preference (converse relation): 
#'  A tuple t1 is better than t2 w.r.t. \code{-p1} if t2 is better than t1 w.r.t. \code{p1}. 
#'  The unary minus operator, i.e. \code{-p1}, is a short hand notation for \code{reverse(p1)}.}
#' }
#' 
#' The function \code{is.complex_pref} returns \code{TRUE} if \code{x} is a complex preference object 
#' (i.e., was constructed by one of these binary operators or the unary operator \code{reverse}) 
#' and \code{FALSE} otherwise.
#' 
#'
#' @seealso See \code{\link{base_pref}} for the construction of base preferences. 
#' See \code{\link{general_pref}} for functions applicable to all kind of preferences.
#' See \code{\link{psel}} for the evaluation of preferences. 
#' 
#' @references 
#' 
#' S. Borzsonyi, D. Kossmann, K. Stocker (2001): The Skyline Operator. In Data Engineering (ICDE '01), pages 421-430.
#' 
#' W. Kiessling (2002): Foundations of Preferences in Database Systems. In Very Large Data Bases (VLDB '02), pages 311-322.
#' 
#' S. Mandl, O. Kozachuk, M. Endres, W. Kiessling (2015): Preference Analytics in EXASolution. 
#' 16th Conference on Database Systems for Business, Technology, and Web.
#' 
#' @examples
#' # define preference for cars with low consumption (high mpg-value) 
#' # and simultaneously high horsepower
#' p1 <- high(mpg) * high(hp)  
#' 
#' # perform the preference search
#' psel(mtcars, p1)
NULL


# Complex preference Construtors
# ------------------------------

# Score functions are NOT needed for C++ algorithms but for igraph, etc.

# Infix Pareto-Constructor
#' @rdname complex_pref
#' @export
"*.preference" <- function(p1, p2) {
  check_pref(p1, p2)
  if (check_empty(p1, p2)) return(get_empty(p1, p2))
  return(paretopref(p1, p2))
}

# Infix Prioritization-Constructor (special constructor as we have to consider prior-chains!)
#' @rdname complex_pref
#' @export
"&.preference" <- function(p1, p2) {
  check_pref(p1, p2)
  if (check_empty(p1, p2)) return(get_empty(p1, p2))
  return(priorpref(p1, p2))
}


# Infix Intersection-Constructor
#' @rdname complex_pref
#' @export
"|.preference" <- function(p1, p2) {
  check_pref(p1, p2)
  if (check_empty(p1, p2)) return(get_empty(p1, p2))
  return(intersectionpref(p1, p2))
}

# Infix-Disjount-Union-Constructor
#' @rdname complex_pref
#' @export
"+.preference" <- function(p1, p2) {
  check_pref(p1, p2)
  if (check_empty(p1, p2)) return(get_empty(p1, p2))
  return(unionpref(p1, p2))
}


# Reverse preference
#' @rdname complex_pref
#' @export
reverse <- function(p) {
  check_pref(p)
  if (is.empty_pref(p)) return(p)
  return(reversepref(p))
}

# This entry has no @rdname as "-" is just used unary!
# (it is exported, but invisible in the documentation
#' @export
"-.preference" <- function(p1, p2) {
  if (nargs() == 1) return(reverse(p1))
  else stop("Operation not defined.")
}


#' @rdname complex_pref
#' @export
is.complex_pref <- function(x) {
  return(inherits(x, "complexpref") || inherits(x, "reversepref"))
}
 
# Helper functions
# ----------------
  
# Helper to check if all given arguments are preferences
check_pref <- function(p1, p2) {
  if (!is.actual.preference(p1) || (nargs() == 2 && !is.actual.preference(p2)))
    stop("This operator requires preference objects as input.")
}

# Check if one (or perhaps both) preference is empty
check_empty <- function(p1, p2) (is.empty_pref(p1) || is.empty_pref(p2))

# Get the result of a complex operation with an empty pref
get_empty <- function(p1, p2) {
  if (is.empty_pref(p1))
    return(p2) # perhaps empty
  else
    return(p1)
}    
