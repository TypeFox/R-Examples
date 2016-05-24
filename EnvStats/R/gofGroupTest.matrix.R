gofGroupTest.matrix <-
function (object, ...) 
{
    nc <- ncol(object)
    if (nc == 1 || !is.numeric(object)) {
        stop(paste("When \"object\" is a matrix, it must", "have 2 or more columns, and it must be numeric."))
    }
    col.names <- colnames(object)
    if (is.null(col.names)) 
        col.names <- paste("Column", 1:nc, sep = ".")
    nr <- nrow(object)
    group <- rep(col.names, each = nr)
    arg.list <- list(object = as.vector(unlist(object)), group = group)
    dots.list <- list(...)
    names.list <- list(data.name = deparse(substitute(object)), 
        group.name = "Columns")
    match.vec <- pmatch(names(dots.list), c("data.name", "group.name"), 
        nomatch = 0)
    if (length(match.vec) == 0 || all(match.vec == 0)) 
        arg.list <- c(arg.list, names.list, dots.list)
    else arg.list <- c(arg.list, list(data.name = deparse(substitute(object)), 
        group.name = "Columns")[-match.vec], dots.list)
    do.call("gofGroupTest.default", args = arg.list)
}
