summaryFull.formula <-
function (object, data = NULL, subset, na.action = na.pass, ...) 
{
    if (missing(object) || (length(object) != 3L)) 
        stop("formula missing or incorrect")
    m <- match.call(expand.dots = FALSE)
    if (is.matrix(eval(m$data, parent.frame()))) 
        m$data <- as.data.frame(data)
    m$... <- NULL
    m$formula <- m$object
    m$object <- NULL
    m$na.action <- na.action
    requireNamespace("stats", quietly = TRUE)
    m[[1L]] <- as.name("model.frame")
    mf <- eval(m, parent.frame())
    response <- attr(attr(mf, "terms"), "response")
    arg.list <- list(object = mf[, response])
    dot.list <- list(...)
    match.vec <- pmatch(names(dot.list), "data.name")
    if (length(match.vec) == 0 || all(is.na(match.vec))) 
        arg.list <- c(arg.list, list(data.name = names(mf)[response]))
    if (ncol(mf) == 1) 
        arg.list <- c(arg.list, dot.list)
    else arg.list <- c(arg.list, list(group = mf[, -response]), 
        dot.list)
    do.call(summaryFull.default, arg.list)
}
