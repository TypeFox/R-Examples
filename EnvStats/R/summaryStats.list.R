summaryStats.list <-
function (object, ...) 
{
    data.name <- deparse(substitute(object))
    length.object <- length(object)
    names.object <- names(object)
    if (is.null(names.object)) 
        names.object <- paste(data.name, "Component", 1:length.object, 
            sep = ".")
    if (all(sapply(object, is.numeric))) {
        if (length.object == 1) {
            arg.list <- list(object = as.vector(unlist(object)))
            match.vec <- pmatch(names(list(...)), "data.name")
            if (length(match.vec) == 0 || is.na(match.vec)) 
                arg.list <- c(arg.list, list(data.name = data.name), 
                  ...)
            else arg.list <- c(arg.list, ...)
            do.call("summaryStats.default", arg.list)
        }
        else {
            group.Ns <- sapply(object, length)
            group <- factor(rep(names.object, group.Ns), levels = names.object)
            summaryStats.default(object = as.vector(unlist(object)), 
                group = group, ...)
        }
    }
    else if (all(sapply(object, is.factor))) {
        list.levels <- lapply(object, levels)
        list.lengths <- sapply(list.levels, length)
        if (!all(list.lengths == list.lengths[1])) 
            stop(paste("When \"object\" is a list and all components are factors,", 
                "all components have to have the same levels"))
        all.levels.identical <- all(sapply(list.levels, function(x) identical(x, 
            list.levels[[1]])))
        if (!all.levels.identical) 
            stop(paste("When \"object\" is a list and all components are factors,", 
                "all components have to have the same levels"))
        group.Ns <- sapply(object, length)
        object <- unlist(object)
        group <- factor(rep(names.object, group.Ns), levels = names.object)
        summaryStats.factor(object = object, group = group, ...)
    }
    else stop(paste("When \"object\" is a list,", "all components must be numeric or all components must be factors"))
}
