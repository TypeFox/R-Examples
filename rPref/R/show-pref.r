#' Show Preferences in Database Query Languages
#' 
#' For a given preference this shows the \code{PREFERRING} clause of a database query in different SQL dialects which support preferences.
#' 
#' @param p A preference.
#' @param dialect The preference query dialect, which determines the syntax of the returned query. 
#'        This has to be one of the following (not case sensitive):
#' \describe{
#'    \item{\code{'EXASOL'}:}{Syntax of the "Skyline" feature of the commercial database EXASOL EXASolution 5.}
#'    \item{\code{'Preference SQL'} or \code{'PSQL'}:}{Syntax of the Preference SQL system. 
#'      This is a research prototype developed at the Chair of Databases and Information Systems of the University of Augsburg. 
#'      See references for details.}
#' }
#' @param df Optional parameter to specify a data frame on which the preference operates causing a partial evaluation. See \code{\link{show.pref}} for details.
#'     
#' @details 
#' 
#' There are few database systems supporting Skyline queries. 
#' A Skyline query consists of a usual SQL query followed by a \code{PREFERRING}-clause (in some rarely used dialects also \code{SKYLINE OF}). 
#' For example consider a database table r(a,b). The preference selection \code{psel(r, low(a) * high(b))} can be expressed by (in the Exasol dialect):
#' 
#' \code{SELECT * FROM r PREFERRING LOW a PLUS HIGH b}
#' 
#' The \code{show.query} function generates just the \code{PREFERRING}-clause, i.e. \code{show.query(low(a) * high(b))} returns 
#' 
#' \code{PREFERRING LOW a PLUS HIGH b}
#' 
#' As usual in SQL queries, all keywords are not case sensitive, i.e., \code{PLUS} or \code{plus} does not make any difference.
#' 
#' @references 
#' 
#' W. Kiessling, M. Endres, F. Wenzel (2011): The Preference SQL System - An Overview. IEEE Data Engineering Bulletin, Vol. 34 No. 3, pages 12-19.  
#' 
#' S. Mandl, O. Kozachuk, M. Endres, W. Kiessling (2015): Preference Analytics in EXASolution. 
#' 16th Conference on Database Systems for Business, Technology, and Web.
#' 
#' @examples
#' 
#' show.query(low(a) * high(b))
#' 
#' show.query(low(a) * high(b), dialect = 'Preference SQL')
#' 
#' @export
show.query <- function(p, dialect = 'EXASOL', df = NULL) {
 
  # Make parameter "dialect" not case sensitive
  dialect <- tolower(dialect)
  
  if (dialect == "psql" || dialect == "preferencesql") dialect <- PSQL
  
  # Check dialect
  if (!(dialect %in% c(EXA, PSQL))) stop(paste0('The dialect "', dialect, '" is unknown.'))
  
  # Check for empty -> return empty string
  if (is.empty_pref(p)) return("")
  
  # Return PREFERRING-clause
  return(paste0('PREFERRING ', show_pref_sql(p, dialect, 'PREF', get_static_terms(df)))) # parent node "PREF" means root node
}


#' Partial Evaluation and String/Expression Output of Preferences
#' 
#' Functions to substitute variables and functions in preferences 
#' which can be calculated before the preference is evaluated on a data frame.
#' This is especially used for the string/expression output of preferences.
#' 
#' @param p The preference to be shown or partially evaluated.
#' @param df (optional) A data frame on which the preference operates.
#' 
#' @details The function \code{pref.str} (or \code{as.character(p)} for a preference \code{p}) returns the preference string 
#' while \code{show.pref} outputs it directly to the console, preceded by \code{'[Preference]'}.
#' If \code{df} is specified, then a partial evaluation of the preference is done before converting it to a string. 
#' With \code{as.expression(p)} the R code for constructing a given preference is given.
#' 
#' The function \code{eval.pref} (with given data frame \code{df}) partially evaluates the internal preference expression and 
#' returns again a preference object. All expressions in \code{p} are evaluated in the environment
#' where \code{p} was defined, except the the column names in \code{df} (which are potential attributes in \code{p}) 
#' and except the special variable \code{df__}, which accesses the entire data set (see \code{\link{psel}}).
#' The content of the data frame \code{df} does not matter; only \code{names(df)} is used to get the "free variables" in \code{p}.
#' 
#' 
#' @section Partial Evaluation Before String Output:
#' 
#' The functions \code{show.pref} and \code{pref.str} have the optional parameter \code{df}.
#' If this parameter is given, these functions call \code{eval.pref} before they output or return the preference string.
#' The following equalities hold:
#' 
#' \itemize{
#'   \item \code{as.character(eval.pref(p, df)) == pref.str(p, df)}
#'   \item \code{show(eval.pref(p, df))} produces the same console output as \code{show.pref(p, df)}
#' }
#' 
#' @seealso See \code{\link{general_pref}} for more utility functions for preferences.
#' 
#' @examples
#' 
#' f <- function(x) 2*x
#' p <- true(cyl == f(1))
#' 
#' # prints 'true(cyl == f(1))'
#' p
#' 
#' # prints 'true(cyl == 2)'
#' show.pref(p, mtcars)
#' eval.pref(p, mtcars)
#' 
#' @export
show.pref <- function(p, df = NULL) {
  # Output string representation on the console
  return(cat('[Preference] ', pref.str(p, df), '\n', sep = ''))
}

# Internal call of evaluation, no validity check of paramters
eval.pref.internal <- function(p, df = NULL) {
  if (is.null(df)) return(p)
  else return(p$evaluate(get_static_terms(df)))
}

#' @export
#' @rdname show.pref
pref.str <- function(p, df = NULL) {
  if (!is.actual.preference(p)) stop("This function needs a preference as first argument.")
  return(eval.pref.internal(p, df)$get_str())
}

#' @export
#' @rdname show.pref
eval.pref <- function(p, df = NULL) {
  if (!is.actual.preference(p)) stop("This function needs a preference as first argument.")
  # Careful: eval(as.expression(p)) would not respect the local environment of the base preferences!
  # We cannot respect the environment of the caller of eval.pref (merging environments did not work here)
  return(eval.pref.internal(p, df))
}

# Get attributes of a data set as symbols
get_static_terms <- function(df) {
  # Get static terms (all attributes of data set, if given)
  if (!is.null(df)) return(lapply(names(df), as.symbol))
  else return(NULL)  
}

# Available Dialects (in lower case)
EXA <- 'exasol'
PSQL <- 'preference sql'

# Translation maps for operators
PSQL_PREF_MAP = list('*'        = 'AND', 
                     '&'        = 'PRIOR TO', 
                     '|'        = 'INTERSECT WITH',
                     '+'        = 'DISJOINT UNION',
                     'reverse'  = 'DUAL')

EXASOL_PREF_MAP = list('*'        = 'PLUS', 
                       '&'        = 'PRIOR TO', 
                       '|'        = NA,
                       '+'        = NA,
                       'reverse'  = 'INVERSE')


# Internal function for converting a preference to its sql representation
show_pref_sql <- function(p, dialect, parent_op = '', static_terms = NULL) {
  
  if (dialect == EXA) use_map <- EXASOL_PREF_MAP
  else                use_map <- PSQL_PREF_MAP
  
  if (is.base_pref(p)) {
    return(show_base_pref(p, dialect, static_terms))
  } else if (is.reversepref(p)) { 
    if (dialect == EXA) # INVERSE is notated as prefix in EXASOL
      return(paste0(use_map[['reverse']], ' (', show_pref_sql(p$p, dialect, static_terms = static_terms), ')'))
    else  # "DUAL" is notated as suffix in Preference SQL!
      return(paste0('(', show_pref_sql(p$p, dialect, static_terms = static_terms), ') ', use_map[['reverse']]))
      
  } else if (is.binarycomplexpref(p)) {
    # usual complex preference (not inverse!)
    opchr <- use_map[[p$op]]
    # Check if operator is available in the given dialect
    if (is.na(opchr)) stop('Operator "', p$op, '" is not available in the dialect "', dialect, '".')
    res <- paste0(show_pref_sql(p$p1, dialect, opchr, static_terms = static_terms), ' ', opchr, ' ', 
                  show_pref_sql(p$p2, dialect, opchr, static_terms = static_terms))
    if (parent_op != opchr && parent_op != 'PREF') res <- paste0('(', res, ')')
    return(res)
  }
}


# Translate Preference to PrefSQL
show_base_pref <- function(p, dialect, static_terms = NULL) {
  
  expr_str <- get_expr_sql(p, static_terms)
  
  if (dialect == EXA) {
    if (is.lowpref(p)) {
      return(paste0('LOW ', expr_str))
    } else if (is.highpref(p)) {
      return(paste0('HIGH ', expr_str))
    } else if (is.truepref(p)) {
      return(paste0(expr_str))
    }
  } else {
    if (is.lowpref(p)) {
      return(paste0(expr_str, ' LOWEST'))
    } else if (is.highpref(p)) {
      return(paste0(expr_str, ' HIGHEST'))
    } else if (is.truepref(p)) {
      return(paste0(expr_str, ' = TRUE'))
    }
  }
}


# Get (evaluated, if static_terms is given) expression of base preference
# Add braces for non-single term 
get_expr_sql <- function(p, static_terms = NULL) {
  # Get expression (eval, if static_terms is given!)
  if (is.null(static_terms)) expr <- p$expr
  else                       expr <- as.expression(p$get_expr_evaled(p$expr, static_terms = static_terms)$expr)
  # Make string of it
  expr_str <- as.character(expr)
  if (length(expr[[1]]) != 1) return(paste0('(', expr_str, ')'))  # embrace complex expressions
  else return(expr_str) 
}
