varGroupTest.formula <-
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
    names.mf <- names(mf)
    data.name <- names.mf[response]
    if (ncol(mf) == 1) {
        group <- rep("Group", length(arg.list$object))
        group.name <- "Group"
    }
    else {
        group <- mf[, -response]
        group.name <- names.mf[-response]
    }
    arg.list <- c(arg.list, list(group = group))
    names.list <- list(data.name = data.name, group.name = group.name)
    dot.list <- list(...)
    match.vec <- pmatch(names(dot.list), c("data.name", "group.name"), 
        nomatch = 0)
    if (length(match.vec) == 0 || all(match.vec == 0)) 
        arg.list <- c(arg.list, names.list, dot.list)
    else arg.list <- c(arg.list, names.list[-match.vec], dot.list)
    if (!missing(data)) 
        arg.list$parent.of.data <- deparse(substitute(data))
    if (!missing(subset)) 
        arg.list$subset.expression <- deparse(substitute(subset))
    do.call(varGroupTest.default, arg.list)
}
