stripChart.formula <-
function (x, data = NULL, dlab = NULL, subset, na.action = NULL, 
    ...) 
{
    if (missing(x) || (length(x) != 3L)) 
        stop("formula missing or incorrect")
    m <- match.call(expand.dots = FALSE)
    if (is.matrix(eval(m$data, parent.frame()))) 
        m$data <- as.data.frame(data)
    m$... <- NULL
    m$formula <- m$x
    m$x <- NULL
    m$na.action <- na.action
    requireNamespace("stats", quietly = TRUE)
    m[[1L]] <- as.name("model.frame")
    mf <- eval(m, parent.frame())
    response <- attr(attr(mf, "terms"), "response")
    if (is.null(dlab)) 
        dlab <- names(mf)[response]
    if (ncol(mf) == 1) 
        arg.list <- list(x = mf[[response]], dlab = dlab, ...)
    else arg.list <- list(x = split(mf[[response]], mf[-response]), 
        dlab = dlab, ...)
    do.call(stripChart.default, arg.list)
}
