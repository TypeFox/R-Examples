gofGroupTest.data.frame <-
function (object, ...) 
{
    nc <- ncol(object)
    if (nc == 1 || !all(sapply(object, is.numeric))) {
        stop(paste("When \"object\" is a data frame,", "it must have 2 or more columns and they must", 
            "all be numeric."))
    }
    names.object <- names(object)
    nr <- nrow(object)
    group <- rep(names.object, each = nr)
    arg.list <- list(object = as.vector(unlist(object)), group = group)
    dots.list <- list(...)
    names.list <- list(data.name = deparse(substitute(object)), 
        group.name = "Columns")
    match.vec <- pmatch(names(dots.list), c("data.name", "group.name"), 
        nomatch = 0)
    if (length(match.vec) == 0 || all(match.vec == 0)) 
        arg.list <- c(arg.list, names.list, dots.list)
    else arg.list <- c(arg.list, names.list[-match.vec], dots.list)
    do.call("gofGroupTest.default", args = arg.list)
}
