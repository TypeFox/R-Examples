
#' Read edits edits from free-form textfile
#'
#' This utility function allows for free editrule definition in a file. One can extract
#' only the numerical (\code{type='num'}), only the categorical (\code{type='cat'}) or all
#' edits (default) in which case an \code{\link{editset}} is returned. 
#' The function first parses all assignments in the file, so it is possible to compute or read
#' a list of categories defining a datamodel for example.
#'
#' @param file name of text file to read in
#' @param type type of edits to extract. Currently, only 'num' (numerical), 'cat' (categorical)  and 'all' are implemented.
#' 
#' @param ... extra parameters that are currently ignored 
#'
#' @return \code{\link{editset}} with all edits if \code{type=all}, \code{\link{editarray}} if \code{type='cat'}, 
#'      \code{\link{editmatrix}} if \code{type='num'}, \code{\link{editset}} with conditional edits if \code{type='mix'}. 
#'   If the return value is a \code{list}, the elements are named \code{numedits} and \code{catedits}.
#'
#' @export
editfile <- function(file,type=c("all","num","cat","mix"), ...){
# TODO: include expandEdits?
    type <- match.arg(type)
    if (!type %in% c('num','cat','all')) stop(paste("type",type,"invalid or not implemented yet"))
    p <- parse(file=file)
    ass <- sapply(p,class) %in% c('<-','=')
    e <- new.env()
    lapply(p[ass],eval,envir=e)
    edits <- p[!ass]
    # substitute constant assignments in rules. The if-statement prevents conversion to list()
    # if the file being read has no actual rules in it.
    if ( length(edits)>0 ){
      edits <- sapply(edits,function(x) as.expression(do.call(substitute,list(x,e))) )
    }
    et <- editTypes(edits)
    numedits <- edits[et == 'num']
    catedits <- edits[et == 'cat']
    mixedits <- edits[et == 'mix']  
    switch(type,
        num = editmatrix(numedits),
        cat = editarray(catedits,env=e),
        mix = editset(mixedits,env=e),
        all = editset(edits,env=e)
    )
}

