calc.relimp.formula <- function (formula, data, weights, na.action, ..., subset = NULL) 
{
### change UG 1.3: added weights argument
    if (missing(formula)) stop("formula missing")
    if (missing(na.action)) 
        na.action <- getOption("na.action")
    m <- match.call(expand.dots = FALSE)
    if (is.matrix(eval(m$data, parent.frame()))) 
        m$data <- as.data.frame(data)
    if (!is.numeric(m$weights)) m$weights <- eval(m$weights, data, parent.frame()) 
    m[[1]] <- as.name("model.frame")
    m$... <- NULL
    mf <- eval(m, parent.frame())
    terms <- attr(mf, "terms")
    resp <- attr(terms, "response")

    if (resp != 1 ) stop("incorrect formula") 
#    if (max(attr(terms, "order")) > 2) 
#        stop("formula contains terms of order higher than 2")
#    if (length(attr(terms,"order")[attr(terms,"order") == 2]) > 2)
#        stop("Up to two 2-factor-interactions are supported only!")
    ## other checks involving other parameters (e.g. type)
    ## done in calc.relimp.lm only
    if (attr(terms, "intercept") != 1) 
        stop("model must contain intercept")
    if (!is.null(attr(mf, "na.action"))) 
        warning(naprint(attr(mf, "na.action")))
    if (!is.null(dim(model.response(mf)))) {
        if (ncol(model.response(mf)) > 1) 
            stop("too many response variables")
    }

    y <- do.call("calc.relimp", list(lm(mf), ...))
    y
}