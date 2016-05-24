#' Write API for the ecdb for a list of basic ecdattr objects
#' 
#' It takes a list of basic ecdattr objects, enrich them in parallel, then save them to ecdb.
#'
#' @method write ecdb
#'
#' @param x a list of basic ecdattr objects
#' @param object an object of ecdb class
#' 
#' @return The row count
#'
#' @keywords ecdb 
#'
#' @export
#'
### <======================================================================>
"write.ecdb" <- function(x, object)
{
    x1 <- x[[1]]
    # print(x1)
    # ecdattr.enrich(x1)
    
    data <- parallel::mclapply(x, ecdattr.enrich, mc.allow.recursive=FALSE)
    
    # print("save data")
    
    ecdattr(object) <- data
    length(data) # row count
}
### <---------------------------------------------------------------------->
#' @rdname write.ecdb
setGeneric("write", function(x, object) standardGeneric("write"))
#' @rdname write.ecdb
setMethod("write", signature("list", "ecdb"), write.ecdb)
### <---------------------------------------------------------------------->
