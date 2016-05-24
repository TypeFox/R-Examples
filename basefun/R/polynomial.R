
polynomial_basis <- function(var, coef, 
                             ui = NULL, ci = NULL) {

    stopifnot(inherits(var, "numeric_var"))
    varname <- variable.names(var)
    support <- support(var)[[varname]]

    stopifnot(all(coef %in% c(0, 1)))
    object <- do.call("polylist", lapply(1:length(coef), function(i) {
        cf <- coef[1:i]
        cf[1:i < i] <- 0
        return(polynomial(cf))
    }))

    if (is.null(ui)) ui <- Diagonal(length(object))
    if (is.null(ci)) ci <- rep(-Inf, nrow(ui))
    stopifnot(nrow(ui) == length(ci))

    basis <- function(data, deriv = 0L) {

        stopifnot(check(var, data))
        if (is.atomic(data)) {
            x <- data
        } else {
            if (is.null(varname)) varname <- colnames(data)[1]
            x <- data[[varname]]
        }
        dobject <- object
        if (deriv > 0) {
            for (i in 1:deriv)
                dobject <- deriv(dobject)
        }
        X <- sapply(dobject, predict, x)
        if (!is.matrix(X)) X <- matrix(X, nrow = 1)
	cn <- c("(Intercept)", varname)
        if (ncol(X) > 2)
            cn <- c(cn, paste(varname, "^", 2:(ncol(X) - 1), sep = ""))
        colnames(X) <- cn[1:ncol(X)]
        if (deriv < 0) X[] <- 0
        attr(X, "constraint") <- list(ui = ui, ci = ci)
        attr(X, "Assign") <- matrix(varname, ncol = ncol(X))
        X
    }

    attr(basis, "variables") <- var
    attr(basis, "intercept") <- TRUE

    class(basis) <- c("polynomial_basis", "basis", class(basis))
    return(basis)
}

model.matrix.polynomial_basis <- function(object, data, deriv = 0L, ...)
    object(data, deriv = .deriv(variable.names(object), deriv))
