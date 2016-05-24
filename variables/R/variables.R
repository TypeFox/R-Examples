
.var <- function(name, desc = NULL) {
    ret <- list(name = name, desc = desc)
    class(ret) <- "var"
    ret
}

factor_var <- function(name, desc = NULL, levels) {
    ret <- .var(name = name, desc = desc)
    ret$support <- factor(levels, levels = levels, labels = levels)
    class(ret) <- c("factor_var", class(ret))
    ret
}

ordered_var <- function(name, desc = NULL, levels) {
    ret <- factor_var(name = name, desc = desc, levels = levels)
    ret$support <- as.ordered(ret$support)
    class(ret) <- c("ordered_var", class(ret))
    ret
}

numeric_var <- function(name, desc = NULL, unit = NULL,
                        support = c(0.0, 1.0), add = c(0, 0), bounds = NULL) {
    ret <- .var(name = name, desc = desc)
    ret$unit <- unit
    stopifnot(length(support) >= 2L)
    stopifnot(all(is.finite(support)))
    stopifnot(is.integer(support) || is.double(support))
    if (is.integer(support) && length(support) == 2L)
        support <- support[1]:support[2]
    ret$support <- support
    discrete <- is.integer(support) || (length(support) > 2L)
    if (discrete) {
        stopifnot(is.null(bounds))
        class(ret) <- c("discrete_var", "numeric_var", class(ret))
        return(ret)
    }
    if (is.null(bounds))
        bounds <- c(-Inf, Inf)
    stopifnot(bounds[1] <= min(support))
    stopifnot(max(support) <= bounds[2])
    ret$bounds <- bounds
    stopifnot(add[1] <= 0 && add[2] >= 0)
    ret$add <- add
    class(ret) <- c("continuous_var", "numeric_var", class(ret))
    ret
}

c.var <- function(...) {
    ret <- list(...)
    nm <- sapply(ret, variable.names)
    stopifnot(all(!duplicated(nm))) ### make sure no duplicate names
    names(ret) <- nm
    stopifnot(all(sapply(ret, function(x) inherits(x, "var"))))
    class(ret) <- "vars"
    ret
}
    
variable.names.var <- function(object, ...)
    object$name

variable.names.vars <- function(object, ...)
    sapply(object, variable.names)

desc <- function(object)
    UseMethod("desc")

desc.var <- function(object)
    object$desc

desc.vars <- function(object)
    sapply(object, desc)

unit <- function(object)
    UseMethod("unit")

unit.numeric_var <- function(object)
    object$unit

unit.var <- function(object)
    return(NA)

unit.vars <- function(object)
    sapply(object, unit)

support <- function(object)
    UseMethod("support")

support.var <- function(object)
    return(structure(list(object$support), 
                     names = variable.names(object)))

support.vars <- function(object)
   structure(do.call("c", lapply(object, support)), 
             names = variable.names(object))

levels.factor_var <- function(x)
    levels(support(x)[[variable.names(x)]])

levels.discrete_var <- function(x)
    support(x)[[variable.names(x)]]

levels.var <- function(x)
    return(NA)

bounds <- function(object)
    UseMethod("bounds")

bounds.continuous_var <- function(object) 
    structure(list(object$bounds),
              names = variable.names(object))

bounds.discrete_var <- function(object) {
    s <- support(object)[[variable.names(object)]]
    structure(list(range(s)), names = variable.names(object))
}

bounds.ordered_var <- function(object) {
    f <- support(object)[[variable.names(object)]]
    structure(list(f[c(1, nlevels(f))]), 
              names = variable.names(object))
}

bounds.vars <- function(object)
   structure(do.call("c", lapply(object, bounds)),
             names = variable.names(object))
    
bounds.default <- function(object)
    structure(list(NA), names = variable.names(object))

is.bounded <- function(object)
    UseMethod("is.bounded")

is.bounded.continuous_var <- function(object)
    any(is.finite(bounds(object)[[variable.names(object)]]))

is.bounded.var <- function(object)
    return(TRUE)

is.bounded.vars <- function(object)
    sapply(object, is.bounded)

mkgrid <- function(object, ...)
    UseMethod("mkgrid")

mkgrid.var <- function(object, ...)
    return(support(object))

mkgrid.continuous_var <- function(object, n = 2, ...) {
    s <- support(object)[[variable.names(object)]]
    add <- object$add
    if (any(max(abs(add)) > 0))
        s <- s + add
    b <- bounds(object)[[variable.names(object)]]
    if (is.finite(b[1]) & (add[1] == 0)) s[1] <- b[1]
    if (is.finite(b[2]) & (add[2] == 0)) s[2] <- b[2]
    stopifnot(n > 0)
    if (n == 1L) return(structure(list(diff(s)), 
                                  names = variable.names(object)))    
    return(structure(list(seq(from = s[1], to = s[2], length.out = n)), 
                     names = variable.names(object))) 
}
    
mkgrid.vars <- function(object, ...)
   structure(do.call("c", lapply(object, mkgrid, ...)),
             names = variable.names(object))

as.data.frame.vars <- function(x, row.names = NULL, optional = FALSE, 
                               n = 1L, ...) {
    g <- mkgrid(x, n = n)
    len <- max(sapply(g, length))
    as.data.frame(lapply(g, function(x) rep_len(x, length.out = len)))
}

as.data.frame.var <- as.data.frame.vars

as.vars <- function(object) 
    UseMethod("as.vars")

as.vars.data.frame <- function(object) {
    v <- lapply(colnames(object), function(x) {
        if (is.ordered(object[[x]])) 
            return(ordered_var(x, levels = levels(object[[x]])))
        if (is.factor(object[[x]])) 
            return(factor_var(x, levels = levels(object[[x]])))
        b <- NULL
        if (is.integer(object[[x]])) {
            s <- sort(unique(object[[x]]))
        } else if (inherits(object[[x]], "Surv")) { ### <FIXME>: only right censored </FIXME>
            s <- c(.Machine$double.eps, max(object[[x]][,1], na.rm = TRUE))
            b <- c(0, Inf)
        } else {
            s <- range(object[[x]], na.rm = TRUE)
        }
        return(numeric_var(x, support = s, bounds = b))
    })
    return(do.call("c", v))
}

check <- function(object, data)
    UseMethod("check")

check.ordered_var <- function(object, data) {
    if (!is.atomic(data)) {
        v <- variable.names(object)
        stopifnot(v %in% names(data))
        data <- data[[v]]
    }
    is.ordered(data) && isTRUE(all.equal(levels(data), 
                                         levels(object)))
}

check.factor_var <- function(object, data) {
    if (!is.atomic(data)) {
        v <- variable.names(object)
        stopifnot(v %in% names(data))
        data <- data[[v]]
    }
    is.factor(data) && isTRUE(all.equal(levels(data), 
                                        levels(object)))
}

check.discrete_var <- function(object, data) {
    if (!is.atomic(data)) {
        v <- variable.names(object)
        stopifnot(v %in% names(data))
        data <- data[[v]]
    }
    all(data %in% support(object)[[variable.names(object)]])
}

check.continuous_var <- function(object, data) {
    if (!is.atomic(data)) {
        v <- variable.names(object)
        stopifnot(v %in% names(data))
        data <- data[[v]]
    }
    b <- bounds(object)[[variable.names(object)]]
    min(data, na.rm = TRUE) >= b[1] && 
    max(data, na.rm = TRUE) <= b[2]
}

check.vars <- function(object, data)
    all(sapply(object, check, data = data))
