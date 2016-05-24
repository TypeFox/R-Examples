weibull.frailty <-
function (formula = formula(data), data = parent.frame(), id = "id", subset, na.action, init, 
        control = list()) {
    call <- match.call()
    m <- match.call(expand.dots = FALSE)
    temp <- c("", "formula", "data", "subset", "na.action")
    m <- m[match(temp, names(m), nomatch = 0)]
    Terms <- if (missing(data)) terms(formula) else terms(formula, data = data)
    m$formula <- Terms
    m[[1]] <- as.name("model.frame")
    m <- eval(m, parent.frame())
    if (NROW(m) == 0) 
        stop("No (non-missing) observations.\n")
    Y <- model.extract(m, "response")
    if (!inherits(Y, "Surv")) 
        stop("Response must be a survival object.\n")
    logT <- log(Y[, 1])
    d <- Y[, 2]
    id <- if (is.character(id) && length(id) == 1) {
        if (missing(data) || !id %in% names(data))
            stop("'id' not a 'data'.\n")
        nam.id <- id
        dd <- if (missing(na.action)) na.omit(data) else na.action(data)
        as.vector(unclass(factor(dd[[id]])))
    } else {
        as.vector(unclass(factor(id)))
    }
    attr(Terms, "intercept") <- 1
    X <- model.matrix(Terms, m)[, -1, drop = FALSE]
    type <- attr(Y, "type")
    if (type != "right") 
        stop("weibull.frailty() currently supports only right-censored data.\n")
    if (missing(init)) 
        init <- NULL
    out <- weibull.frailty.fit(logT, d, X, id, init, control)
    out$y <- Y
    out$x <- X
    out$id <- id
    out$nam.id <- nam.id
    out$terms <- Terms
    out$data <- if (missing(data)) m else data
    out$call <- call
    class(out) <- "weibull.frailty"
    out
}
