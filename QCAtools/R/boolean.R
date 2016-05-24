#' And
#'
#' Logical 'and' of two conditions
#' @param v1 A vector of fuzzy set scores of cases
#' @param v2 A vector of fuzzy set scores of cases
#' @return the fuzzy set scores of the logical conjunction of v1 and v2 for each case, i.e.
#' the minimum in each component
#' 
#' @rdname boolean
#' 
#' @examples
#' and(c(0,0.5,1), c(0.25, 0.75, 0.75))
#' or(c(0,0.5,1), c(0.25, 0.75, 0.75))
#' not(c(0,0.5,1))
#' 
#' @export
and <- function(v1, v2) {
  return(pmin(v1, v2))
}

#' Or
#'
#' Logical 'or' of two conditions
#' @return the fuzzy set scores of the logical disjunction of v1 and v2 for each case, i.e.
#' the maximum in each component
#' @rdname boolean
#' @export
or <- function(v1, v2) {
  return(pmax(v1, v2))
}

#' Not
#'
#' Logical 'not' of a condition
#' @param v A vector of fuzzy set scores of cases
#' @return the fuzzy set scores of the negation of v for each case, i.e. 1-v
#' @rdname boolean
#' @export
not <- function(v) {
  return(1 - v)
}

#' Convert formula to function
#'
#' When given a Boolean formula (in disjunctive normal form, see details), this
#' function produces a function that takes a data.frame of a QCA data
#' table and computes the fuzzy set score for each case of membership
#' in the set described by the formula
#' 
#' @details
#' Set names must be capitalized in the formula and the data; if they are
#' lowercase, they are interpreted as the negation of the set.
#' If \code{formula} is a string, logical 'or' is expressed as a '+',
#' and logical 'and' as a '*'.
#' If \code{formula} is a list of strings, the strings are assumed to be
#' the dosjuncts and are concatenated with '+'.
#' Disjunctive normal form means that the formula must be a disjunction of
#' conjunctions of elementary or negated elementary sets. Example:
#' \code{A*b*C + a*B}
#' 
#' 
#' @param formula A string or vector of strings containing a Boolean formula in disjunctive normal form
#' @return a function that takes a data.frame and computes the fuzzy set score
#' of the set described by the formula for each case into a vector
#' 
#' @examples
#' formula_to_function("A*b*C + a*B")
#' 
#' @export
formula_to_function <- function(formula) {
  if (length(formula) > 1) {
    paths <- formula
  } else {
    paths <- stringr::str_split(formula, "\\s*\\+\\s*")[[1]]    
  }
  conjuncts <- vector(length=length(paths))
  for (i in 1:length(paths)) {
    conjuncts[i] <- list(stringr::str_split(paths[i], "\\s*\\*\\s*")[[1]])
  }
  return(function(data) {
    names(data) <- toupper(names(data))
    results <- rep(0, times=nrow(data))
    for (i in 1:length(paths)) {
      conjresults <- rep(1, times=nrow(data))
      for(j in 1:length(conjuncts[[i]])) {
        
        condname <- toupper(conjuncts[[i]][[j]])
        newcond <- data[,condname]
        if (condname != conjuncts[[i]][[j]]) { ## lower case name
          newcond <- not(newcond)
        }
        conjresults <- and(newcond, conjresults)
      }
      results <- or(results, conjresults)
    }
    return(results)
  })
}

#' Evaluate a formula
#'
#' When given a Boolean formula (see details) and a \code{data.frame} of cases and fuzzy
#' set score for conditions, computes for each case the score of the membership
#' in the set described by the formula
#' 
#' @details
#' If \code{formula} is a function, it must take a \code{data.frame} and return
#' a vector.
#' 
#' If \code{formula} is a string or list of strings, the following conventions hold:
#' Set names must be capitalized in the formula and the data; if they are
#' lowercase, they are interpreted as the negation of the set.
#' If \code{formula} is a string, logical 'or' is expressed as a '+',
#' and logical 'and' as a '*'.
#' If \code{formula} is a list of strings, the strings are assumed to be
#' the dosjuncts and are concatenated with '+'.
#' The formula must be in disjunctive normal form, i.e. it must be a disjunction of
#' conjunctions of elementary or negated elementary sets. Example:
#' \code{A*b*C + a*B}
#' 
#' 
#' @param formula A string, list of strings or function representing a Boolean formula in disjunctive normal form
#' @param data A data frame where the rows represent cases and the columns the sets. Column names must be as in the formula.
#' @return the fuzzy set score of the set described by the formula for each case in the data
#' 
#' @examples
#' require(QCAGUI)
#' data(d.urban)
#' evaluate_dnf(d.urban, "MLC*frb + CP")
#' 
#' @export
evaluate_dnf <- function(data, formula) {
  if (!is.function(formula)) {
    formula <- formula_to_function(formula)
  }
  return(formula(data))
}

#' Compute the consistency value
#'
#' Compute a consistency score for an implication/necessity/sufficiency statement.
#' 
#' @description
#' Computes the consistency score of "formula1 -> formula2" (sufficient condition) 
#' or "formula1 <- formula2" (necessary condition), depending on whether \code{type}
#' is "->" or "<-".
#' If \code{type} is "<->" it computes an equivalence score of formula1 and formula2
#' via the formula \code{sum(min(X,Y))/(sum(max(X,Y))}
#' 
#' @details
#' If \code{formula} is a function, it must take a \code{data.frame} and return
#' a vector.
#' 
#' If \code{formula} is a string or list of strings, the following conventions hold:
#' Set names must be capitalized in the formula and the data; if they are
#' lowercase, they are interpreted as the negation of the set.
#' If \code{formula} is a string, logical 'or' is expressed as a '+',
#' and logical 'and' as a '*'.
#' If \code{formula} is a list of strings, the strings are assumed to be
#' the dosjuncts and are concatenated with '+'.
#' The formula must be in disjunctive normal form, i.e. it must be a disjunction of
#' conjunctions of elementary or negated elementary sets. Example:
#' \code{A*b*C + a*B}
#' 
#' 
#' @param formula1 A string, list of strings or function representing a Boolean formula in disjunctive normal form
#' @param formula2 A string, list of strings or function representing a Boolean formula in disjunctive normal form
#' @param type either "->", "<-" or "<->", depending on the direction of the implication that is to be evaluated
#' @param data A data frame where the rows represent cases and the columns the sets. Column names must be as in the formula.
#' @return the consistency score of the implication described by \code{formula1}, \code{type} and \code{formula2}
#' 
#' @examples
#' require(QCAGUI)
#' data(d.urban)
#' consistency("MLC + FRB", "->", "CP", d.urban)
#' 
#' @export
consistency <- function(formula1, type="->", formula2, data) {
  col1 <- evaluate_dnf(data, formula1)
  col2 <- evaluate_dnf(data, formula2)
  if (type == "->") {
    return(sum(pmin(col1,col2))/sum(col1))
  } else if (type == "<-") {
    return(sum(pmin(col1,col2))/sum(col2))
  } else if (type == "<->") {
    return(sum(pmin(col1,col2))/sum(pmax(col1,col2)))
  }
}
