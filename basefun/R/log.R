
log_basis <- function(var, ui = c("none", "increasing", "decreasing"), 
                      remove_intercept = FALSE) {

    stopifnot(inherits(var, "numeric_var"))
    varname <- variable.names(var)
    support <- support(var)[[varname]]
    stopifnot(support[1] > 0)
    bounds <- bounds(var)[[varname]]
    if (is.finite(bounds[1]))
        stopifnot(bounds[1] >= 0)

    ui <- match.arg(ui)
    ci <- switch(ui, "none" = -Inf, 0)
    ui <- switch(ui, "none" = matrix(1),
                     "increasing" = matrix(1),
                     "decreasing" = matrix(-1))
    if (!remove_intercept) {
        ci <- c(-Inf, ci)
        ui <- cbind(0, ci)
    }

    basis <- function(data, deriv = 0L) {

        stopifnot(check(var, data))
        if (is.atomic(data)) {
            x <- data
        } else {
            if (is.null(varname)) varname <- colnames(data)[1]
            x <- data[[varname]]
        }
        ### x <= 0 not allowed
        x <- pmax(x, .Machine$double.eps)
        if (deriv == 0) {
            X <- cbind(1, matrix(log(x), ncol = 1))
        } else if (deriv == 1) {
            X <- cbind(0, matrix(1 / x, ncol = 1))
        } else if (deriv > 1) {
            stop("not yet implemented")
        } else if (deriv < 0) {
            X <- cbind(0, matrix(0, ncol = 1, nrow = length(x)))
        }
        colnames(X) <- c("(Intercept)", paste("log(", varname, ")", sep = ""))
        if (remove_intercept) {
            X <- X[, -1L, drop = FALSE]
        } else {
            varname <- c("(Intercept)", varname)
        }
        attr(X, "constraint") <- list(ui = ui, ci = ci)
        attr(X, "Assign") <- varname
        X
    }

    attr(basis, "variables") <- var
    attr(basis, "intercept") <- !remove_intercept

    class(basis) <- c("log_basis", "basis", class(basis))
    return(basis)
}

model.matrix.log_basis <- function(object, data, deriv = 0L, ...)
    object(data, deriv = .deriv(variable.names(object), deriv))
