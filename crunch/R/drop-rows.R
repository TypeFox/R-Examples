##' Permanently delete rows from a dataset
##'
##' @param dataset a \code{CrunchDataset}
##' @param expr a \code{CrunchLogicalExpr}
##' @return \code{dataset} without the rows indicated by \code{expr}
##' @seealso \code{\link{exclusion}} for a non-destructive way to suppress rows
##' @examples
##' \dontrun{
##' ds <- dropRows(ds, ds$gender == "Male")
##'}
##' @export
dropRows <- function (dataset, expr) {
    ## Construct a ZCL delete payload, then POST it
    payload <- list(command="delete",
        filter=zcl(expr))
    out <- crPOST(shojiURL(dataset, "fragments", "table"), body=toJSON(payload))
    dropCache(self(dataset))
    return(dataset)
}
