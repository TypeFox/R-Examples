kendallTrendTest.formula <-
function (y, data = NULL, subset, na.action = na.pass, ...) 
{
    if (missing(y) || (length(y) != 3L)) 
        stop("formula missing or incorrect")
    m <- match.call(expand.dots = FALSE)
    if (is.matrix(eval(m$data, parent.frame()))) 
        m$data <- as.data.frame(data)
    m$... <- NULL
    m$formula <- m$y
    m$y <- NULL
    m$na.action <- na.action
    requireNamespace("stats", quietly = TRUE)
    m[[1L]] <- as.name("model.frame")
    mf <- eval(m, parent.frame())
    response <- attr(attr(mf, "terms"), "response")
    arg.list <- list(y = mf[, response])
    names.mf <- names(mf)
    names.list <- list(data.name = names.mf[response])
    if (ncol(mf) > 1) {
        x <- mf[, -response]
        arg.list <- c(arg.list, list(x = x))
        names.list <- c(names.list, list(data.name.x = names.mf[-response]))
    }
    dot.list <- list(...)
    match.vec <- pmatch(names(dot.list), c("data.name", "data.name.x"), 
        nomatch = 0)
    if (length(match.vec) == 0 || all(match.vec == 0)) 
        arg.list <- c(arg.list, names.list, dot.list)
    else arg.list <- c(arg.list, names.list[-match.vec], dot.list)
    if (!missing(data)) 
        arg.list$parent.of.data <- deparse(substitute(data))
    if (!missing(subset)) 
        arg.list$subset.expression <- deparse(substitute(subset))
    do.call(kendallTrendTest.default, arg.list)
}
