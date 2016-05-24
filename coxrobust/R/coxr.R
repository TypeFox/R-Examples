
coxr <- function(formula, data, subset, na.action, trunc = 0.95,
                f.weight = c("linear", "quadratic", "exponential"),
                singular.ok = TRUE, model = FALSE) {

    call <- match.call()

    mf <- match.call(expand.dots = FALSE)
    m  <- match(c("formula", "data", "subset", "na.action"), names(mf),
                nomatch = 0)
    mf <- mf[c(1, m)]
    mf[[1]] <- as.name("model.frame")
    mf <- eval(mf, parent.frame())
    if ( NROW(mf) == 0 ) {
        stop("no (non-missing) observations")
    }

    mterms <- attr(mf, "terms")

    y <- model.extract(mf, "response")
    if ( !inherits(y, "Surv") ) {
        stop("response must be a \"Surv\" object")
    } else {
        type <- attr(y, "type")
        if ( type != "right" ) {
           stop(sprintf("\"%s\" type of survival data is not supported", type))
        }
    }
    
    x <- model.matrix(mterms, mf)
    x <- x[, -1, drop = FALSE]

    if ( missing(trunc) ) {
       trunc <- 0.95
    }

    if ( missing(f.weight) ) {
        f.weight <- "quadratic"
    } else {
        f.weight <- match.arg(f.weight)
    }
    f.weight = which(f.weight == c("linear", "quadratic", "exponential"))

    init <- rep(0,ncol(x))
    fit <- coxr.fit(x, y, trunc, init, f.weight, singular.ok)

    if ( length(fit$skip) > 0 ) {

        skipcol <- paste(fit$skip, collapse = ", ")

        msg <- sprintf("X matrix deemed to be singular; variable %s",
                        skipcol)

        if ( !singular.ok || length(fit$skip) == ncol(x) ) {
            stop(msg)
        } else {
            warning(msg)
        }

    }
    
    class(fit) <- "coxr"

    fit$na.action <- attr(mf, "na.action")
    fit$call <- call
    fit$terms <- mterms

    if (model) {
        fit$model <- mf
    }

    fit$x <- x
    fit$y <- y

    fit

}