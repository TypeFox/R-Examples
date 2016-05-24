`Stratiplot.formula` <- function(formula,  data, subset,
                                 na.action = "na.pass",
                                 type = "l", ylab = NULL, xlab = "",
                                 pages = 1, ...) {
    cl <- match.call()
    mf <- match.call(expand.dots = FALSE)
    m <- match(c("formula", "data", "subset", "na.action"),
               names(mf), 0L)
    mf <- mf[c(1L, m)]
    mf$drop.unused.levels <- TRUE
    mf[[1L]] <- as.name("model.frame")
    ## force eval of the na.action argument default
    mf$na.action <- substitute(na.action)
    mf <- eval(mf, parent.frame())
    mt <- attr(mf, "terms")
    y <- model.response(mf, "numeric")
    data <- data.frame(model.matrix(mt, mf), check.names = FALSE)[,-1]
    names(data) <- gsub("`", "", names(data))
    n.vars <- ncol(data)
    y <- rep(y, n.vars)
    if(is.null(ylab))
        ylab <- as.character(attr(mt, "variables")[[2]])
    Stratiplot.default(x = data, y = y, type = type,
                       ylab = ylab, xlab = xlab, pages = pages,
                       ...)
}
