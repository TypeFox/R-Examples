
weights.mlt <- function(object, ...) {
    if (!is.null(object$weights))
        return(object$weights)
    rep(1, NROW(object$data))
}

coef.mlt <- function(object, fixed = TRUE, ...) {
    if (fixed) 
        return(object$coef)
    return(object$par)
}

"coef<-" <- function(object, value)
    UseMethod("coef<-")

"coef<-.mlt" <- function(object, value) {
    cf <- coef(object, fixed = TRUE)
    stopifnot(length(cf) == length(value))
    if (!is.null(names(value)))
        stopifnot(all.equal(names(cf), names(value)))
    object$par <- object$parm(value)
    object$coef[] <- value ### preserve names
    object
}

### use maxLik::hessian?
Hessian <- function(object, ...)
    UseMethod("Hessian")

Hessian.mlt <- function(object, parm = coef(object, fixed = FALSE), ...)
    object$hessian(parm, weights = weights(object))
    
Gradient <- function(object, ...)
    UseMethod("Gradient")

Gradient.mlt <- function(object, parm = coef(object, fixed = FALSE), ...)
    as.vector(colSums(estfun(object, parm = parm)))

vcov.mlt <- function(object, parm = coef(object, fixed = FALSE), ...)
    solve(Hessian(object, parm = parm))

logLik.mlt <- function(object, parm = coef(object, fixed = FALSE), 
                       w = weights(object), ...) {
    ret <- -object$loglik(parm, weights = w)
    ###    attr(ret, "df") <- length(coef(object, fixed = FALSE))
    attr(ret, "df") <- object$df
    class(ret) <- "logLik"
    ret
}

estfun.mlt <- function(object, parm = coef(object, fixed = FALSE), 
                       w = weights(object), ...)
    -object$score(parm, weights = w)

mkgrid.mlt <- function(object, n, ...)
    mkgrid(object$model, n = n, ...)

mkgrid.ctm <- function(object, n = n, ...)
    mkgrid(object$model, n = n, ...)

variable.names.mlt <- function(object, ...)
    variable.names(object$model)

model <- function(object)
    UseMethod("model")

.one_factor_only <- function(object) {
    f <- inherits(object, "formula_basis")
    v <- as.vars(object)
    f && (length(v) == 1 && inherits(v[[1L]], "factor_var"))
}
    
model.ctm <- function(object) {
    x <- object$bases
    ret <- list(response_trafo = c("continuous", "discrete")[.one_factor_only(x$response) + 1L],
                response_type = sapply(x$response, class),
                response_var = as.vars(x$response),
                interaction_trafo = !is.null(x$interacting),
                shift_trafo = !is.null(x$shifting))
    if (ret$interaction_trafo) {
        ret$interaction_vars = as.vars(x$interacting)
        ret$interaction_type = c("continuous", "discrete")[.one_factor_only(x$interacting) + 1L]
    }
    if (ret$shift_trafo) 
        ret$shift_vars = as.vars(x$shifting) 
    return(ret)
}

model.mlt <- function(object) {
    c(model(object$model), 
      list(todistr = object$todistr$name,
           fixed = object$fixed))
}

description <- function(object) {
    stopifnot(inherits(object, "mlt"))
    m <- model(object)
    cond <- m$interaction_trafo || m$shift_trafo
    strat <- m$interaction_trafo && m$interaction_type == "discrete"
    lin <- cond && (!m$interaction_trafo || strat)
    if (lin) pm <- switch(m$todistr, "logistic" = "odds",
                                     "minimum extreme value" = "hazards", "")

    ret <- paste(m$response_trafo, 
       if (!cond) "unconditional",
       if (strat) "stratified",
       if (lin) "linear", 
       if (!lin && cond) "conditional",
       "transformation model", 
       "(transformed", m$todistr, "distribution)")
    if (cond && lin && pm != "")
        ret <- c(ret, paste(m$response_trafo, if (strat) "stratified", "proportional", pm, "model"))
    ret <- gsub("\\s\\s*", " ", ret)
    return(ret)
}

summary.mlt <- function(object, ...) {

    ret <- list(call = object$call,
                convergence = object$convergence,
                type = paste(description(object), collapse = "\n\t"),
                logLik = logLik(object),
                AIC = AIC(object),
                coef = coef(object))
    class(ret) <- "summary.mlt"
    ret
}

print.summary.mlt <- function(x, digits = max(3L, getOption("digits") - 3L), ...) {

    cat("\nCall:\n")
    print(x$call)
    if (x$convergence != 0L)
    cat("\nCould not estimate parameters; optimisation did not converge!\n")
    cat("\nType: ", x$type)
    cat("\nAIC: ", x$AIC)
    cat("\nLog-Likelihood: ", x$logLik, " (df = ", attr(x$logLik, "df"), ")", sep = "")
    cat("\n")
    cat("\nCoefficients:", x$coef)
    cat("\n\n")
    invisible(x)
}

print.mlt <- function(x, ...)
    print(summary(x, ...))

as.vars.ctm <- function(object)
    as.vars(object$model)

as.vars.mlt <- function(object)
    as.vars(object$model)

bounds.ctm <- function(object)
    bounds(as.vars(object))

bounds.mlt <- function(object)
    bounds(as.vars(object))

print.response <- function(x, ...) {

    ac <- as.character
    obs <- paste(ifelse(!is.na(x$exact), ac(x$exact), 
                 paste("(", ac(x$cleft), ", ", ac(x$cright), "]", sep = "")))

    if (all(is.na(x$tleft) & is.na(x$tright))) {
        print(obs, quote = FALSE, ...)
        return(invisible(obs))
    }

    trc <- character(length(obs))
    i <- (!is.na(x$tleft) & is.na(x$tright))
    if (sum(i) > 0)
        trc[i] <- paste("| >", x$tleft[i])
    i <- (is.na(x$tleft) & !is.na(x$tright))
    if (sum(i) > 0)
        trc[i] <- paste("| <", x$tright[i])
    i <- (!is.na(x$tleft) & !is.na(x$tright))
    if (sum(i) > 0)
        trc[i] <- paste("| (", x$tleft[i], ", ", x$tright[i], "]", sep = "")
    ret <- paste("{", obs, trc, "}", sep = "")
    print(ret, quote = FALSE, ...)
}
