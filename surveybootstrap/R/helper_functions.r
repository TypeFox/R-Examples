#######################################################################
## helper_functions.R
##
## various helper functions for use in the networkreporting
## package
##

##########################################################################
##' get a variable from a dataframe or vector
##'
##' this function was written because a few of the estimator functions
##' need to use weights, and there are several cases to handle:
##' the user could pass in a column name, a vector of weights, or
##' nothing (in which case, the weights should default to 1 for each
##' row in the dataset). for the special case of getting weights, look
##' at the curried fn get.weights (right below)
##' 
##' @param survey.data the survey dataset
##' @param var either NULL, a column name, or a vector of values
##' @param default the default value to fill in if the variable
##'        is not found
##' @return a vector of values whose length is the same as the
##'         number of rows in survey.data; if var is NULL, this has
##'         the default values
##' @keywords internal
get.var <- function(survey.data, var, default=NA) {

  ## weights will default to 1 for everyone, unless the user specified
  ## a weights variable
  if (is.null(var)) {
    
    return(rep(default, nrow(survey.data)))
    
  } else if (length(var) == 1) {
  
    ## ... otherwise, see if the weights variable is referring
    ## to a column of the dataframe; try to
    ## grab sampling weights from survey dataframe
    var.vals <- try(subset(survey.data,
                           select=var),
                    silent=TRUE)

    if( inherits(var.vals, "try-error") ||
       ncol(var.vals) != 1 ||
       ! is.numeric(var.vals[,1]) ) {

      stop(paste(var,
                 " does not identify a valid column in the data.\n"))
    }

    var <- as.numeric(var.vals[,1])

    ## NB: this is a workaround -- may want to rethink this in the future.
    ## see
    ## http://stackoverflow.com/questions/21618423/extract-a-dplyr-tbl-column-as-a-vector
    if (inherits(var, "tbl_df")) {
      var <- collect(select(var, 1))[[1]]
    }

    return(var)
    
  } else if (is.numeric(var) & length(var) == nrow(survey.data)) {

    ## if var a vector with one entry per row, then these
    ## are our values
    return(var)
  } else {    
    stop("can't determine what the values should be for ", var, ".")
  }
    
  
}

##########################################################################
##' get the weights column from a dataframe
##'
##' this is the same as get.var with the default value set to 1
##' instead of NA
##' @param ... (this is a function curried from \code{get.var})
##' @keywords internal
get.weights <- functional::Curry(get.var, default=1)

##########################################################################
##' grab a function based on its name
##'
##' helper to grab a fn that is passed in as an argument
##' 
##' this is based on Hadley Wickham's response to an SO
##' post: \url{http://stackoverflow.com/questions/14183766/match-fun-provide-error-with-functions-defined-inside-functions}
##' with some minor modifications
##'
##' @param fn the function to search for
##' @param env the environment to start searching in
##' @return fn, if fn is already a function; otherwise, the first function found
##'         in env or one of its parents whose name is fn
##' @keywords internal 
get.fn <- function(fn, env = parent.frame()) {

    ## base case: fn is already a function
    if (is.function(fn)) {
      return(fn)
    }
  
    ## base case: nothing left to search through
    if (identical(env, emptyenv())) {
        stop("Could not find function ", fn, "!")
    }

    ## base case: found function in env
    if (exists(fn, env, inherits=FALSE) &&
        is.function(env[[fn]])) {
        return(env[[fn]])

    ## recursive case: look through the environment
    ## above env
    } else {
        return(get.fn(fn, parent.env(env)))
    }

}

##########################################################################
##' get a variable from a dataframe or vector
##'
##' this function was written because a few of the estimator functions
##' need to use weights, and there are several cases to handle:
##' the user could pass in a column name, a vector of weights, or
##' nothing (in which case, the weights should default to 1 for each
##' row in the dataset). for the special case of getting weights, look
##' at the curried fn get.weights (right below)
##'
##' @param survey.data the survey dataset
##' @param var either NULL, a column name, or a vector of values
##' @param default the default value to fill in if the variable
##'        is not found
##' @return a vector of values whose length is the same as the
##'         number of rows in survey.data; if var is NULL, this has
##'         the default values
##' @keywords internal
get.var <- function(survey.data, var, default=NA) {

  ## weights will default to 1 for everyone, unless the user specified
  ## a weights variable
  if (is.null(var)) {

    return(rep(default, nrow(survey.data)))

  } else if (length(var) == 1) {

    ## ... otherwise, see if the weights variable is referring
    ## to a column of the dataframe; try to
    ## grab sampling weights from survey dataframe
    var.vals <- try(survey.data[,var,drop=FALSE],
                    silent=TRUE)

    ## if( inherits(var.vals, "try-error") ||
    ##    ncol(var.vals) != 1 ||
    ##    ! is.numeric(var.vals[,1]) ) {
    if( inherits(var.vals, "try-error") ||
       ncol(var.vals) != 1) {

      stop(paste(var,
                 " does not identify a valid column in the data.\n"))
    }

    var <- var.vals[,1]

    if (inherits(var.vals, "tbl_df")) {
        var <- collect(select(var.vals, 1))[[1]]
    }

    return(var)

  } else if (length(var) == nrow(survey.data)) {

    ## if var a vector with one entry per row, then these
    ## are our values
    return(var)
  } else {
    stop("can't determine what the values should be for ", var, ".")
  }


}

##########################################################################
##' get the weights column from a dataframe
##'
##' this is the same as get.var with the default value set to 1
##' instead of NA
##' @param ... (this is a function curried from \code{get.var})
##' @keywords internal
get.weights <- functional::Curry(get.var, default=1)

##########################################################################
##' only prints things out in verbose mode
##'
##' @param verbose if TRUE, print things out; otherwise, do nothing
##' @param ... arguments to pass to cat if verbose is TRUE
##' @keywords internal
vcat <- function(verbose=TRUE, ...) {

  if(verbose) {
    message(...)
  }

  invisible()
}

##########################################################################
##' parse a formula that describes the design of a survey
##'
##' Given a formula of the form\cr
##' \code{~ psu_v1 + psu_v2 + ... + strata(strata_v1 + strata_v2 + ...)}\cr
##' \itemize{
##'  \item{"psu.formula"}{a formula describing the primary sampling unit vars}
##'  \item{"strata.formula"}{a formula describing the strata (if any)}
##' }\cr
##' TODO
##' \itemize{
##'  \item{}{check to be sure no response is included (or warn)}
##'  \item{}{check formulas for strata more carefully...}
##' }
##'
##' @param formula a formula describing the sample design (see above)
##' @return a list with entries \code{psu.formula} and \code{strata.formula}
##' @keywords internal
parse_design <- function(formula) {

  ## see http://stackoverflow.com/questions/10224805/how-to-select-a-part-of-formula-in-formula-in-r
  ## for some helpful info

  psu.formula <- formula
  strata.formula <- NULL

  these.labels <- attr(terms(formula), "term.labels")

  strata.idx <- grep("strata\\(", these.labels)

  if (length(strata.idx) == 1) {

    # grab the expression in the strata(...) part of the formula
    #strata.text <- str_match

    strata.text <- stringr::str_match(these.labels[strata.idx],
                                      "strata\\((.+)\\)")[2]

    ## updating instead of creating a new formula b/c this preserves
    ## the environment that the original formula was created in...
    strata.formula <- update(formula,
                             paste("~ ", strata.text))

    psu.formula <- update.formula(formula,
                                  paste("~ . - strata(",strata.text,")"))

  } else if (length(strata.idx > 1)) {

    stop("Cannot have more than one strata() specification in the design formula.")
  }

  return(list(psu.formula=psu.formula,
              strata.formula=strata.formula))

}


##########################################################################
##' compute the weighted mean
##'
##' given a vector of values and a vector of weights, compute the
##' weighted mean
##'
##' @param x the vector of values
##' @param w the vector of weights
##' @param na.rm if TRUE, only consider elmeents of x that are not missing
##'              (and their corresponding entries in w). Defaults to FALSE.
##' @return the weighted mean
##' @keywords internal
weighted.mean <- function(x, w, na.rm=FALSE) {

  if (na.rm) {
    idx <- (1:length(x))[!is.na(x)]
  } else {
    idx <- 1:length(x)
  }

  return(sum(x[idx]*w[idx])/sum(w[idx]))
}

