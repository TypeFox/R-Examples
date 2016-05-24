#' Loglinear Model Utilities

#' These functions generate lists of terms to specify a loglinear model
#' in a form compatible with loglin and provide for conversion to an
#' equivalent loglm specification.  They allow for a more conceptual
#' way to specify such models.


#' models of joint independence, of some factors wrt one or more other factors

#' @param nf number of factors for which to generate model
#' @param table a contingency table used for factor names, typically the output from \code{\link[base]{table}}
#' @param factors names of factors used in the model when \code{table} is not specified
#' @param with    indices of the factors against which others are considered jointly independent
#' @export 

joint <- function(nf, table=NULL, factors=1:nf, with=nf) {
    if (!is.null(table)) factors <- names(dimnames(table))
    if (nf == 1) return (list(term1=factors[1]))
    if (nf == 2) return (list(term1=factors[1], term2=factors[2]))
    others <- setdiff(1:nf, with)
    result <- list(term1=factors[others], term2=factors[with])
    result
  }


#' models of conditional independence of some factors wrt one or more other factors

#' @param nf number of factors for which to generate model
#' @param table a contingency table used for factor names, typically the output from \code{\link[base]{table}}
#' @param factors names of factors used in the model when \code{table} is not specified
#' @param with    indices of the factors against which others are considered conditionally independent
#' @export 

conditional <- function(nf, table=NULL, factors=1:nf, with=nf) {
    if (!is.null(table)) factors <- names(dimnames(table))
    if (nf == 1) return (list(term1=factors[1]))
    if (nf == 2) return (list(term1=factors[1], term2=factors[2]))
    main <- setdiff(1:nf, with)
    others <- matrix(factors[with], length(with), length(main))
    result <- rbind(factors[main], others)
    result <- as.list(as.data.frame(result, stringsAsFactors=FALSE))
    names(result) <- paste('term', 1:length(result), sep='')
    result
  }

#' models of mutual independence of all factors

#' @param nf number of factors for which to generate model
#' @param table a contingency table used for factor names, typically the output from \code{\link[base]{table}}
#' @param factors names of factors used in the model when \code{table} is not specified
#' @export 

mutual <- function(nf, table=NULL, factors=1:nf) {
   if (!is.null(table)) factors <- names(dimnames(table))
   result <- sapply(factors[1:nf], list)
   names(result) <- paste('term', 1:length(result), sep='')
   result  
  }


#' saturated model: highest-order interaction

#' @param nf number of factors for which to generate model
#' @param table a contingency table used for factor names, typically the output from \code{\link[base]{table}}
#' @param factors names of factors used in the model when \code{table} is not specified
#' @export 

saturated <- function(nf, table=NULL, factors=1:nf) {
	if (!is.null(table)) factors <- names(dimnames(table))
	list(term1=factors[1:nf])
}

# models of conditional independence, given one pair of variables
## Not needed: handled by condit, with length(with)>1

#condit2 <- function(nf, factors=1:nf, with=1:2) {
#    if (nf == 1) return (list(term1=factors[1]))
#    if (nf == 2) return (list(term1=factors[1], term2=factors[2]))
#    others <- setdiff(1:nf, with)
#    result <- rbind(factors[with], cbind(factors[others], factors[others]))
#    result <- as.list(as.data.frame(result, stringsAsFactors=FALSE))
#    names(result) <- paste('term', 1:length(result), sep='')
#    result
#}

#' markov models of a given order

#' @param nf number of factors for which to generate model
#' @param table a contingency table used for factor names, typically the output from \code{\link[base]{table}}
#' @param factors names of factors used in the model when \code{table} is not specified
#' @param order   order of the markov chain
#' @export 

markov <- function(nf, factors=1:nf, order=1) {
    if (nf == 1) return (list(term1=factors[1]))
    if (nf == 2) return (list(term1=factors[1], term2=factors[2]))
    if (length(factors) < order+2) {
      warning(paste('Not enough factors for order', order, 'Markov chain; using order=1'))
      order <-1
      result <- rbind(factors[1:(nf-1)], factors[2:nf])
    }
    else {
      if (nf <= order+1) result <- factors[1:nf]
      else {
        result <- NULL
        for (i in 1:(order+1)) 
          result <- rbind(result, factors[i:(nf-order+i-1)])
      }
    }

    result <- as.list(as.data.frame(result, stringsAsFactors=FALSE))
    names(result) <- paste('term', 1:length(result), sep='')
    result
}

#' convert a loglin model to a model formula for loglm

#' @param  x a list of terms in a loglinear model, such as returned by \code{joint}, \code{conditional}, \dots
#' @param  env environment in which to evaluate the formula
#' @source Code from Henrique Dallazuanna, <wwwhsd@gmail.com>, R-help 7-4-2013

loglin2formula <- function(x, env = parent.frame()) {
	terms <- lapply(x, paste, collapse = ":")
	formula(sprintf(" ~ %s", do.call(paste, c(terms, sep = "+"))), env=env)
}

#' convert a loglin model to a string, using bracket notation for the high-order terms

#' @param x a list of terms in a loglinear model, such as returned by \code{joint}, \code{conditional}, \dots
#' @param brackets characters to use to surround model terms.  Either a single character string containing two characters
#'        or a character vector of length two.
#' @param sep characters used to separate factor names within a term
#' @param collapse  characters used to separate terms
#' @param abbrev

loglin2string <- function(x, brackets = c('[', ']'), sep=',', collapse=' ', abbrev) {
	if (length(brackets)==1 && (nchar(brackets)>1)) brackets <- unlist(strsplit(brackets, ""))
	terms <- lapply(x, paste, collapse=sep)
	terms <- paste(brackets[1], terms, brackets[2], sep='')
	paste(terms, collapse= ' ')  
}


