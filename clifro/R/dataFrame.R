# dataFrame Class ---------------------------------------------------------

#' @importFrom methods setClass
setClass("dataFrame", slots = c(names = "character",
                                row.names = "character"),
         contains = "list")

#' @rdname dim
#' @aliases dimnames,dataFrame-method
#' @importFrom methods setMethod
setMethod("dimnames", "dataFrame",
          function (x) 
          {
            list(x@row.names, names(x))
          }
)

#' Dimension Attributes of a Clifro Object
#' 
#' Retrieve the dimensions or dimension names of a \code{dataFrame} object.
#' 
#' @seealso \code{\link{cf_query}} for creating \code{cfData} objects, and
#'   \code{\link{cf_station}} for creating \code{cfStation} objects.
#' 
#' @param x a \code{dataFrame} object
#' 
#' Specifically, a \code{dataFrame} object is any \code{\link{cfStation}} or 
#' \code{cfData} object. These functions are provided for the user to have (some)
#' familiar \code{data.frame}-type functions available for use on \pkg{clifro}
#' objects.
#' 
#' @rdname dim
#' @aliases dim,dataFrame-method
#' @importFrom methods setMethod
setMethod("dim", "dataFrame",
          function (x) 
          {
            c(length(x@row.names), length(x))
          }
)

#' @importFrom methods setMethod
setMethod("row.names", "dataFrame",
          function (x) 
          {
            x@row.names
          }
)

#' @importFrom methods setMethod
setMethod("rownames", "dataFrame",
          function (x, do.NULL = TRUE, prefix = "row") 
          {
            dn <- dimnames(x)
            if (!is.null(dn[[1L]])) 
              dn[[1L]]
            else {
              nr <- NROW(x)
              if (do.NULL) 
                NULL
              else if (nr > 0L) 
                paste0(prefix, seq_len(nr))
              else character()
            }
          }
)

#' @importFrom methods setMethod
setMethod("colnames", "dataFrame",
          function (x, do.NULL = TRUE, prefix = "col") 
          {
            dn <- dimnames(x)
            if (!is.null(dn[[2L]])) 
              dn[[2L]]
            else {
              nc <- NCOL(x)
              if (do.NULL) 
                NULL
              else if (nc > 0L) 
                paste0(prefix, seq_len(nc))
              else character()
            }
          }
)

#' @importFrom methods as setMethod
#' @rdname Extract
#' @aliases [[,dataFrame-method
setMethod("[[", "dataFrame",
          function (x, i) 
            as(x, "data.frame")[[i]]
)

#' @importFrom methods as setMethod
#' @rdname Extract
#' @aliases [,dataFrame,ANY,ANY,ANY-method
setMethod("[", "dataFrame", 
          function (x, i, j, drop)
            as(x, "data.frame")[i, j, drop]
)

#' @importFrom methods setMethod
#' @rdname Extract
#' @aliases $,dataFrame-method
setMethod("$", "dataFrame",
          function (x, name)
          {
            which_col = pmatch(name, names(x))
            if (!is.na(which_col)){
              if (is.na(match(name, names(x))))
                warning("name partially matched in dataFrame")
              return(x[, which_col])
            }
            NULL
          }
)

#' @importFrom methods as setMethod
setMethod("show", "dataFrame",
          function (object) 
          {
            print(as(object, "data.frame"))
          }
)

#' @importFrom stats setNames
#' @importFrom methods setAs
setAs("dataFrame", "data.frame",
      function(from){
        setNames(data.frame(from@.Data), from@names)
      })

#' @importFrom methods as setMethod
setMethod("as.data.frame", "dataFrame",
          function (x)
          {
            as(x, "data.frame")
          }
)