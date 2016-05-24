#' UrlData -- unified access to WWW resources
#' 
#' This class provides the infrastructure to
#' scrape the web with a Extract, Transform, Load (ETL)
#' approach.
#'
#' In most cases, it is not necessary to subclass \code{UrlData}.
#' The slots can be set by the \code{urldata} function and allow
#' to customize each step of the process.
#'
#' @seealso \code{\link{urldata}}
#' 
#' @examples
#' getSlots("UrlData")
#' 
#' @exportClass UrlData
#' @name UrlData-class
#' @rdname UrlData-class
setClass(
  Class="UrlData", 
  representation=representation(
    template="character", 
    resource="character",
    map.lst="list", 
    extract.fct="function", 
    transform.fct="function"
  ),
  contains="Xdata"
)

#' Constructor for UrlData objects
#'
#' This function creates a web resource. As this, it
#' allows to map an Web API into an R S4 class.
#'
#' @param resource        the name of the resource. Required.
#' @param template        a pattern for the url. Must contain \%s for substitution. Required.
#' @param extract.fct     a function that takes an URI and returns the raw data. Default readLines.
#' @param transform.fct   a function that takes the raw data and returns the cleaned/transformed data. Default identity.
#' @param clss            name of the class to create. Default UrlData, must be inherited from this class.
#' @param ...             parameters for the query. Must be named arguments, values can be characters (for defaults), 
#'                        NULL, or functions.
#'
#' @export
urldata <- function(resource, template, extract.fct=readLines, transform.fct=identity, clss="UrlData", ...) {
    new(clss, resource=resource, template=template, map.lst=list(...), extract.fct=extract.fct, transform.fct=transform.fct)
}

#' @rdname query-methods
#' @name query
#' @export
#' @docType methods
#' @aliases query query,UrlData,character-method
setMethod(
  f="query",
  signature=c(self="UrlData", resource="character"),
  definition=function(self, resource, verbose=getOption("verbose"), ...) {
    if(resource==self@resource) {
        if(verbose) cat("constructing URL..\n")
        mapped <- list()
        arg.lst <- list(...)
        unused <- setdiff(names(arg.lst), names(self@map.lst))
        if(length(unused)>0) warning("unused argument(s) to query: '", paste(unused, collapse="', '"), "'")
        for(n in names(self@map.lst)) {
            arg.fct <- self@map.lst[[n]]
            arg <- arg.lst[[n]]
            if(!is.null(arg)) 
                mapped[[n]] <- if(is.function(arg.fct)) arg.fct(arg) else as.character(arg)
            else
                mapped[[n]] <- if(is.function(arg.fct)) arg.fct() else arg.fct
        }
        na.arg <- sapply(mapped, is.na)
        if(any(na.arg)) stop("query with missing required parameters: '", paste(names(which(na.arg)), collapse="', '"), "'")
        uri <- strsubst(self@template, mapped)
    
        if(verbose) cat("Downloading ", uri,"...\n")
        res <- try(self@extract.fct(uri), silent=TRUE)
        if(inherits(res, "try-error")) return(res)
    
        if(verbose) cat("Transforming data ...\n")
        res <- try(self@transform.fct(res), silent=TRUE)
        return(res)
    }
    
    if(verbose) cat("trying inherited method..\n")
    callNextMethod(self=self, resource=resource, verbose=verbose, ...)
  }
)

#' @rdname queries-methods
#' @name queries
#' @export
#' @docType methods
#' @aliases queries queries,UrlData-method
setMethod(
  f="queries",
  signature="UrlData",
  definition=function(self) { 
    ret <- c(callNextMethod(), self@resource)
    names(ret) <- NULL
    ret <- ret[ret!="character"]
    return(ret)
  }
)
