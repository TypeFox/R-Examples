summaryStats.data.frame <-
function (object, ...) 
{
    if (all(sapply(object, is.numeric))) {
        names.object <- names(object)
        nc <- ncol(object)
        if (nc == 1) {
            arg.list <- list(object = as.vector(unlist(object)))
            match.vec <- pmatch(names(list(...)), "data.name")
            if (length(match.vec) == 0 || is.na(match.vec)) 
                arg.list <- c(arg.list, list(data.name = names.object), 
                  ...)
            else arg.list <- c(arg.list, ...)
            do.call("summaryStats.default", arg.list)
        }
        else {
            nr <- nrow(object)
            group <- rep(names.object, each = nr)
            summaryStats.default(object = as.vector(unlist(object)), 
                group = group, ...)
        }
    }
    else if (all(sapply(object, is.factor))) {
        list.levels <- lapply(object, levels)
        list.lengths <- sapply(list.levels, length)
        if (!all(list.lengths == list.lengths[1])) 
            stop(paste("When \"object\" is a data frame and all columns are factors,", 
                "all columns have to have the same levels"))
        all.levels.identical <- all(sapply(list.levels, function(x) identical(x, 
            list.levels[[1]])))
        if (!all.levels.identical) 
            stop(paste("When \"object\" is a data frame and all columns are factors,", 
                "all columns have to have the same levels"))
        names.object <- names(object)
        nr <- nrow(object)
        object <- unlist(object)
        group <- rep(names.object, each = nr)
        summaryStats.factor(object = object, group = group, ...)
    }
    else stop(paste("When \"object\" is a data frame,", "all columns must be numeric or all columns must be factors"))
}
