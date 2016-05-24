
as.basis.formula <- function(object, data = NULL, remove_intercept = FALSE, 
                             ui = NULL, ci = NULL, negative = FALSE, scale = FALSE,
                             ...) {

    if (inherits(data, "data.frame")) {
        vars <- as.vars(data[all.vars(object)])
    } else {
        stopifnot(!scale)
        if (inherits(data, "var")) data <- c(data)
        vars <- do.call("c", data[all.vars(object)])
        data <- as.data.frame(vars, n = 10)
    }
    varnames <- all.vars(object)
    stopifnot(varnames %in% variable.names(vars))

    mf <- NULL
    if (!is.null(data) && NROW(data) > 0) {
        ### see http://developer.r-project.org/model-fitting-functions.txt
        mf <- model.frame(object, data = data, drop.unused.levels = TRUE)
        mt <- attr(mf, "terms")
        X <- model.matrix(mt, data = mf, ...)
        contr <- attr(X, "contrasts")
        xlevels <- .getXlevels(mt, mf)
        if (scale) {
            mins <- apply(X, 2, min, na.rm = TRUE)
            maxs <- apply(X, 2, max, na.rm = TRUE)
            if (attr(mt, "intercept") == 1)
                mins["(Intercept)"] <- 0
            maxs <- maxs - mins
        }
    }

    ret <- function(data, deriv = 0L) {

        if (!is.null(vars)) ### intercept only
            stopifnot(check(vars, data))
        data <- data[varnames]
        stopifnot(is.data.frame(data)) 
        if (!is.null(mf)) {
            mf <- model.frame(mt, data = data, xlev = xlevels)
            if(!is.null(cl <- attr(mt, "dataClasses"))) .checkMFClasses(cl, mf)
            X <- model.matrix(mt, mf, contrasts = contr)
        } else {
            mf <- model.frame(object, data)
            X <- model.matrix(attr(mf, "terms"), data = mf, ...)
        }
        if (scale) {
            X <- X - matrix(mins, nrow = nrow(X), ncol = length(mins), byrow = TRUE)
            X <- X / matrix(maxs, nrow = nrow(X), ncol = length(mins), byrow = TRUE)
        }
        if (remove_intercept) {
            a <- attr(X, "assign")
            X <- X[, colnames(X) != "(Intercept)", drop = FALSE]
            attr(X, "assign") <- a[-1]
        }
        if (is.null(ui)) ui <- Diagonal(ncol(X))
        if (is.null(ci)) ci <- rep(-Inf, nrow(ui))
        stopifnot(nrow(ui) == length(ci))
        attr(X, "constraint") <- list(ui = ui, ci = ci)
        attr(X, "Assign") <- c("(Intercept)", varnames)[attr(X, "assign") + 1]
        if (deriv == 1) {
            nm <- names(deriv)[deriv == 1L]
            X[, attr(X, "Assign") %in% nm] <- 1
            X[, !(attr(X, "Assign") %in% nm)] <- 0
        }
        if (deriv < 0)
            X[] <- 0
        if (negative) return(-X)
        return(X)
    }

    attr(ret, "variables") <- vars
    ### note: ~ a - 1 also contains intercept!
    attr(ret, "intercept") <- !remove_intercept 
    class(ret) <- c("formula_basis", "basis", class(ret))
    ret
}

model.matrix.formula_basis <- function(object, data, dim = NULL, deriv = 0L, ...) {

    if (!is.null(dim)) {
        nd <- names(dim)
        nd <- nd[nd %in% variable.names(object)]
        if (length(nd) > 1 & all(dim[nd] > 1)) {
            data <- do.call("expand.grid", data[nd])
        } else {
            if (length(nd) > 1 & (sum(dim[nd] > 1) > 1))
                stop("either all or just one element of dim can be larger one")
            stopifnot(length(unique(sapply(data[variable.names(object)], 
                                           length))) == 1)
            data <- as.data.frame(data[variable.names(object)])
        }
    }

    object(data, deriv = .deriv(variable.names(object), deriv))
}
