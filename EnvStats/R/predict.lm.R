predict.lm <-
function (object, ...) 
{
    ret.val <- stats::predict.lm(object, ...)
    dots.list <- list(...)
    if (length(dots.list)) {
        names.dots.list <- names(dots.list)
        index <- !is.na(pmatch(names.dots.list, "se.fit"))
        if (any(index)) {
            se.fit <- as.logical(dots.list[[names.dots.list[index]]])
            if (se.fit) 
                ret.val <- c(ret.val, n.coefs = object$rank)
        }
    }
    ret.val
}
