
Legendre_basis <- function(var, order = 2, 
                           ui = c("none", "increasing", "decreasing", "cyclic",
                                  "positive", "negative"), ...) { # "zerointegral")

    stopifnot(inherits(var, "numeric_var"))

    object <- legendre.polynomials(order, ...)

    ui <- match.arg(ui)
    B <- Bernstein_basis(var, order, ui = ui)
    constr <- get("constr", environment(B))

    varname <- variable.names(var)
    support <- range(support(var)[[varname]])
    stopifnot(all(diff(support) > 0))

    basis <- function(data, deriv = 0L, integrate = FALSE) {

        stopifnot(check(var, data))
        if (integrate) stop("Integration not yet implemented")

        stopifnot(order > max(c(0, deriv)))
        if (is.atomic(data)) {
            x <- data
        } else {
            if (is.null(varname)) varname <- colnames(data)[1]
            x <- data[[varname]]
        }
        dobject <- object
        if (deriv > 0) {
            for (i in 1:deriv)
                dobject <- sapply(dobject, deriv)
        }

        ### map into [0, 1]
        x <- (x - support[1]) / diff(support)
        stopifnot(all(x >= 0 && x <= 1))

        ### Legendre is on [-1, 1] !
        X <- do.call("cbind", as.vector(lapply(dobject, predict, 2 * x - 1)))
        colnames(X) <- paste("L", 1:ncol(X), sep = "")
        if (deriv > 0)
            X <- X * (2 / diff(support)^deriv)
        if (deriv < 0) X[] <- 0
        if (ui != "none")
            constr$ui <- constr$ui %*% L2B(order)
        attr(X, "constraint") <- constr
        attr(X, "Assign") <- matrix(varname, ncol = ncol(X))
        X
    }

    attr(basis, "variables") <- var
    attr(basis, "intercept") <- TRUE

    class(basis) <- c("Legendre_basis", "Bernstein_basis", "basis", class(basis))
    return(basis)
}
