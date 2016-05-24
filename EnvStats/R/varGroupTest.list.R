varGroupTest.list <-
function (object, ...) 
{
    length.object <- length(object)
    if (length.object == 1 || !all(sapply(object, is.numeric))) 
        stop(paste("When 'object' is a list,", "it must have 2 or more components, and", 
            "all components must be numeric"))
    names.object <- names(object)
    if (is.null(names.object)) 
        names.object <- paste("Component", 1:length.object, sep = ".")
    group.Ns <- sapply(object, length)
    group <- factor(rep(names.object, group.Ns), levels = names.object)
    arg.list <- list(object = as.vector(unlist(object)), group = group)
    dots.list <- list(...)
    names.list <- list(data.name = deparse(substitute(object)), 
        group.name = "Components")
    match.vec <- pmatch(names(dots.list), c("data.name", "group.name"), 
        nomatch = 0)
    if (length(match.vec) == 0 || all(match.vec == 0)) 
        arg.list <- c(arg.list, names.list, dots.list)
    else arg.list <- c(arg.list, list(data.name = deparse(substitute(object)), 
        group.name = "Columns")[-match.vec], dots.list)
    do.call("varGroupTest.default", args = arg.list)
}
