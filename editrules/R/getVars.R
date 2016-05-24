
#' get names of variables in a set of edits
#'
#' @param E \code{\link{editset}}, \code{\link{editmatrix}}, or \code{\link{editarray}}
#' @param ... Arguments to be passed to or from other methods
#' @seealso \code{\link{getA}}, \code{\link{getb}}, \code{\link{getAb}}, \code{\link{getOps}}
#' @example ../examples/getVars.R
#' @return \code{character} vector with the names of the variables. 
#' @export
getVars <- function(E,...){
    UseMethod("getVars")
}


#' Returns the variable names of an (in)equality \code{editmatrix} E
#'
#' @export
#' @method getVars editmatrix
#' @keywords internal
getVars.editmatrix <- function(E,...){
  colnames(E)[-ncol(E)]
}

#' Returns the variable names of an (in)equality \code{editmatrix} E
#'
#' @param type should unique variable names, colnames, all variable names or category names be extracted?
#' @export
#' @method getVars cateditmatrix
#' @keywords internal
getVars.cateditmatrix <- function(E, type=c("uniquevar", "colnames","var", "cat"), ...){
  nms <- colnames(E)[-ncol(E)]
  var <- sub(":.+", "", nms)
  cat <- sub(".+:", "", nms)
  cat[var==cat] <- "TRUE"
  switch( match.arg(type)
        , colnames = nms
        , var = var
        , cat = cat
        , unique(var)
        )
}

#' get variable names in editarray
#'
#' @export
#' @method getVars editarray
#' @keywords internal
getVars.editarray <- function(E,type='cat',...){ 
   if (!type=='cat') return(NULL)
   names(attr(E,"ind"))
}

#' getr variable names
#'
#' @method getVars editset
#' @param type (editset- or list only) select which variables to return. \code{all} means all (except dummies), \code{num} means 
#'      all numericals, \code{cat} means all categoricals, \code{mix} means those numericals appearing in a logical 
#'      constraint and \code{dummy} means dummy variables connecting the logical with numerical constraints.
#' @export
#' @rdname getVars
getVars.editset <- function(E, type=c('all','num','cat','mix','dummy'), ...){
    type <- match.arg(type)
    numvars <- c()
    catvars <- c()

    if (type %in% c('all','num')){
        numvars <- unique(c(getVars(E$num), getVars(E$mixnum)))
    }
    if ( type == 'mix' ) numvars <- getVars(E$mixnum)
    if ( type %in% c('all','cat')){
        catvars <- getVars(E$mixcat)
        catvars <- catvars[!catvars %in% rownames(E$mixnum)]
    }
    if ( type == 'dummy'){
        catvars <- rownames(E$mixnum)
    }
    c(numvars, catvars)
}

#' get variable names
#' @method getVars NULL
#' @export
#' @rdname getVars
getVars.NULL <- function(E,...){
    NULL
}

#' get variable names
#' @export
#' @method getVars editlist
#' @keywords internal
getVars.editlist <- function(E,...){
# under normal circumstances, each part of an editlist has the same variables
    if ( length(E) == 0 ) return(NULL)
    getVars.editset(E[[1]], ...)
}


