summaryFull.data.frame <-
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
            do.call("summaryFull.default", arg.list)
        }
        else {
            nr <- nrow(object)
            group <- rep(names.object, each = nr)
            summaryFull.default(object = as.vector(unlist(object)), 
                group = group, ...)
        }
    }
    else stop(paste("When \"object\" is a data frame,", "all columns must be numeric"))
}
