
##' Sort data according to columns in data frame
##'
##' @title Sort data frame
##' @param data Data frame
##' @param x variable to order by
##' @param ... additional variables to order by
##' @return data.frame
##' @export
##' @examples
##' data(hubble)
##' dsort(hubble, "sigma")
##' dsort(hubble, hubble$sigma,"v")
##' dsort(hubble,~sigma+v)
dsort <- function(data,x,...) {
    if (missing(x)) return(data)
    if (inherits(x,"formula")) x <- all.vars(x)
    if (is.character(x) && length(x)<nrow(data)) x <- lapply(x,function(z) data[,z])
    dots <- list(...)
    args <- lapply(dots, function(x) {
        if (length(x)==1 && is.character(x)) x <- data[,x]
        x
    })
    if (!is.list(x)) x <- list(x)
    data[do.call("order",c(x,args)),]
}
