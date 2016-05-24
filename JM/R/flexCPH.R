flexCPH <-
function (formula = formula(data), data = parent.frame(), subset, na.action, init, control = list()) {
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
    con <- list(lng.in.kn = 3, ord = 4, knots = NULL, numeriDeriv = "fd", eps.Hes = 1e-06, parscale = NULL)
    con[names(control)] <- control
    Y <- model.extract(m, "response")
    if (!inherits(Y, "Surv")) 
        stop("Response must be a survival object.\n")
    logT <- log(Y[, 1])
    d <- Y[, 2]
    attr(Terms, "intercept") <- 1
    X <- model.matrix(Terms, m)
    X <- X[, -1, drop = FALSE]
    type <- attr(Y, "type")
    if (type != "right") 
        stop("flexCPH() supports currently only right-censored data.\n")
    if (missing(init)) 
        init <- NULL
    out <- flexCPH.fit(logT, d, X, init, con)
    out$control <- con
    out$terms <- Terms
    out$call <- call
    class(out) <- "flexCPH"
    out
}
