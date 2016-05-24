#' Create a \code{ffdf} from All Combinations of Factors
#'
#' Similar as \code{expand.grid} in the base package generates an \code{ffdf}.
#' Code is almost copy-pasted from \code{\link[base]{expand.grid}}.
#'
#' @export
#' @example ../examples/expand_ffgrid.R
#' @param ... \code{ff} vectors, \code{ff} factors or a list containing these.
#' @param KEEP.OUT.ATTRS currently ignored
#' @param stringsAsFactors logical specifying if character vectors are converted to factors. Irrelevant for \code{ff} as character vectors are factors in 
#' package ff.
#' @return A \code{ffdf} containing one row for each combination of the supplied factors. The first factors vary fastest. 
#' The columns are labelled by the factors if these are supplied as named arguments or named components of a list. 
#' @seealso \code{\link[base]{expand.grid}}
expand.ffgrid <- function (..., KEEP.OUT.ATTRS = TRUE, stringsAsFactors = TRUE) {
    nargs <- length(args <- list(...))
    if(nargs == 1){
      if(inherits(args[[1]], "list")){
        args <- args[[1]]
        nargs <- length(args)
      }
    }
    #if (!nargs) 
    #    return(as.data.frame(list()))
    #if (nargs == 1L && is.list(a1 <- args[[1L]])) 
    #    nargs <- length(args <- a1)
    #if (nargs == 0L) 
    #    return(as.data.frame(list()))
    cargs <- vector("list", nargs)
    iArgs <- seq_len(nargs)
    nmc <- paste("Var", iArgs, sep = "")
    nm <- names(args)
    if (is.null(nm)) 
        nm <- nmc
    else if (any(ng0 <- nzchar(nm))) 
        nmc[ng0] <- nm[ng0]
    names(cargs) <- nmc
    rep.fac <- 1L
    d <- sapply(args, length)
    #if (KEEP.OUT.ATTRS) {
    #    dn <- vector("list", nargs)
    #    names(dn) <- nmc
    #}
    orep <- prod(d)
    if (orep == 0L) {
        for (i in iArgs) cargs[[i]] <- args[[i]][FALSE]
    }
    else {
        for (i in iArgs) {
            x <- args[[i]]
            #if (KEEP.OUT.ATTRS) 
            #    dn[[i]] <- paste(nmc[i], "=", if (is.numeric(x)) 
            #      format(x)
            #    else x, sep = "")
            nx <- length(x)
            orep <- orep/nx
            x <- x[ffrep.int(ffrep.int(ffseq_len(nx), ffrep.int(rep.fac, nx)), orep)] 
            #x <- x[rep.int(rep.int(seq_len(nx), rep.int(rep.fac, nx)), orep)]
            #if (stringsAsFactors && !is.factor(x) && is.character(x)) 
            #    x <- factor(x, levels = unique(x))
            cargs[[i]] <- x
            rep.fac <- rep.fac * nx
        }
    }
    #if (KEEP.OUT.ATTRS) 
    #    attr(cargs, "out.attrs") <- list(dim = d, dimnames = dn)
    #rn <- .set_row_names(as.integer(prod(d)))
    #structure(cargs, class = "data.frame", row.names = rn)
    as.ffdf.list(cargs)
}
