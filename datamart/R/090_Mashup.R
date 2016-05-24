#' Combine resource into a mashup object
#'
#' The \code{Mashup} class administers a list of Xdata objects. This
#' can be data objects representing different data sources
#' such as internal data or web data. It can also be 
#' calculated data resources, so-called resource functions
#' of class \code{ResFunc}.
#'
#' In a way, this class can be viewed as a make-like tool
#' for data. The resource functions can declare dependencies. 
#' When a resource is requested by the \code{query} method,
#' the \code{Mashup} class takes care of the build order.
#'
#' @examples
#' getSlots("Mashup")
#'
#' @seealso \code{\link{datamart}}
#'
#' @exportClass Mashup
#' @name Mashup-class
#' @rdname datamart
setClass(
  Class="Mashup", 
  representation=representation(res.env="environment"),
  contains="Xdata"
)

#' Constructor for Mashup objects
#'
#' @param ...   named arguments of Xdata objects
#' @param clss  name of the class to create. Default Mashup, must be inherited from this class.
#'
#' @export
#' @rdname datamart
datamart <- function(..., clss="Mashup") {
    # some checks
    arg.lst <- list(...)
    if(length(arg.lst)==0) stop("datamart() needs at least one parameter.")
    
    # populate resource data structure
    res.env <- new.env()
    for (o in arg.lst) {
        if(!inherits(o, "Xdata")) stop("invalid argument to datamart(), needs to derive from Xdata.")
        res.lst <- queries(o)
        
        # no resource?
        if(length(res.lst)==0) stop("data object with no resource passed to datamart(): '", as.character(o), "'.")
        
        for (r in res.lst) {
            # resource already present?
            if(exists(r, envir=res.env, inherits=FALSE)) 
                stop("double resource entry: '", r, "'.")
            
            # result
            assign(r, o, envir=res.env)
        }
    }
    
    # check if all dependencies are present
    r <- ls(envir=res.env)
    r1 <- unlist(lapply(r, function(s) dependencies(res.env[[s]])))
    missed <- setdiff(r1, r)
    if(length(missed)>0) stop("depended resources not found: '", paste(missed, collapse="', '"), "'.", sep="")

    # build object 
    new(clss, res.env=res.env)
}

#' @rdname query-methods
#' @name query
#' @export
#' @docType methods
#' @aliases query query,Mashup,character-method
setMethod(
    f="query",
    signature=c(self="Mashup", resource="character"),
    definition=function(self, resource, verbose=TRUE, ...) {
        # build dependency vector
        check_deps <- resource
        deps <- c()
        while(length(check_deps)>0) {
            r <- check_deps[[1]]
            check_deps <- tail(check_deps, -1)
            newdeps <- unlist(dependencies(self@res.env[[r]]))
            check_deps <- unique(c(newdeps, check_deps))
            deps <- c(newdeps, deps)
        }
        deps <- unique(deps)
        if(verbose) cat("Mashup: build order for resource '", resource, "': '", paste(deps, collapse="', '"), "'.\n", sep="")
        
        # build resources and return the last one
        # TODO: make it optionally persistent
        if(length(deps)>0) {
            cached <- new.env()
            for (r in deps) {
                o <- self@res.env[[r]]
                subdeps <- dependencies(self@res.env[[r]])
                if(verbose) {
                    if(length(subdeps)>0) 
                        cat("Mashup builds '", r, "' and passes pre-built resources '", paste(subdeps, collapse="', '"), "'.\n", sep="")
                    else
                        cat("Mashup builds '", r, "' which requires no pre-built resources.\n", sep="")
                }
                args <- list()
                args[["self"]] <- o
                args[["resource"]] <- r
                args[["verbose"]] <- verbose
                for (d in subdeps) args[[d]] <- cached[[d]]
                ret <- do.call(query, args)
                assign(r, ret, envir=cached)
            }
        }

        # provide the resource requested
        if(!exists(resource, envir=self@res.env)) stop("Don't know how to devliver '",resource,"'.", sep="")
        if(verbose) cat("Mashup builds requested resource '", resource, "'.\n", sep="")
        args <- list(...)
        args[["self"]] <- self@res.env[[resource]]
        args[["resource"]] <- resource
        args[["verbose"]] <- verbose
        for (d in dependencies(self@res.env[[resource]])) args[[d]] <- cached[[d]]
        do.call(query, args)
    }
)

#' @rdname queries-methods
#' @name queries
#' @export
#' @docType methods
#' @aliases queries queries,Mashup-method
setMethod(
  f="queries",
  signature="Mashup",
  definition=function(self) { 
    qrs <- ls(envir=self@res.env)
    ret <- c(callNextMethod(), qrs)
    names(ret) <- NULL
    return(ret)
  }
)
