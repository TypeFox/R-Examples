readDimensionSDML <- function(x)
{
    dim <- NULL
    dimnames <- NULL
    names <- list()
    if (xmlName(x) == "dimension") {
        if (xmlSize(x)) {
            ind <- 0
            for (k in 1:xmlSize(x)) {
                if (xmlName(x[[k]]) != "text") {
                    dim[ind + 1] <- xmlAttrs(x[[k]])["size"]
                    if (!is.na(xmlAttrs(x[[k]])["name"]))
                        dimnames[ind + 1] <- xmlAttrs(x[[k]])["name"]
                    names[[ind + 1]] <- getDataSDML(xmlChildren(x[[k]]))
                    ind <- ind + 1
                }
            }
        }
        if (ind == 0)
            dim <- 0
    }
    mode(dim) <- "integer"
    if (!is.null(dimnames)) names(names) <- dimnames
    list(dim = dim, names = names)
}
