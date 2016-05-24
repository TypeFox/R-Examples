
### box product of multiple _named_ basi(e)s objects
b <- function(..., sumconstr = FALSE) {

    bases <- list(...)
    stopifnot(all(sapply(bases, inherits, what = c("basis", "bases"))))
    bnames <- names(bases)
    stopifnot(length(unique(bnames)) == length(bases))

    varnames <- sapply(bases, variable.names)
    stopifnot(all(!is.null(varnames)))

    attr(bases, "sumconstr") <- sumconstr

    class(bases) <- c("box_bases", "bases")
    bases
}

model.matrix.box_bases <- function(object, data, model.matrix = TRUE,
    dim = NULL, deriv = NULL, integrate = NULL, ...) {

    varnames <- variable.names(object)
    bnames <- names(object)

    if (model.matrix & !is.data.frame(data))
        data <- expand.grid(data[varnames])

    if (!is.null(deriv)) {
        stopifnot(length(deriv) == 1)
        ### if (!names(deriv) %in% varnames) NULL
    }
    if (!is.null(integrate)) {
        stopifnot(length(integrate) == 1)
        if (!names(integrate) %in% varnames) integrate <- NULL
    }
    ret <- lapply(bnames, function(b) {
        thisargs <- list()
        thisargs$object <- object[[b]]
        thisargs$data <- data
        if (!is.null(deriv)) {
            if (names(deriv) %in% variable.names(object[[b]]) || 
                !(names(deriv) %in% varnames))  ### all X need to be matrix(0) then
                thisargs$deriv <- deriv
            ### don't feed deriv to other basis functions
            ### because we look at their _product_ here
        }
        if (!is.null(integrate)) {
            if (names(integrate) %in% variable.names(object[[b]]))
                thisargs$integrate <- integrate
        }
        if (!is.null(dim))
            thisargs$dim <- dim[names(dim) %in% variable.names(object[[b]])]
        X <- do.call("model.matrix", thisargs)
        attr(X, "Assign") <- rbind(attr(X, "Assign"), b)
        X
    })
    if (!model.matrix) return(ret)
    if (attr(object, "sumconstr")) {
        s1 <- mkgrid(object[[1]], 2)
        X1 <- model.matrix(object[[1]], data = expand.grid(s1))
        s2 <- mkgrid(object[[2]], 2) ### min/max for numerics; all levels for factors
        if (any(sapply(s2, function(x) is.integer(x))))
            warning("integer-valued variable with sumcontr = TRUE")
        X2 <- model.matrix(object[[2]], data = expand.grid(s2))
        if (min(X2, na.rm = TRUE) < 0)
            stop("negative values in model matrix not allowed for sumconstr = TRUE")
        ui <- attr(X1, "constraint")$ui
        ci <- attr(X1, "constraint")$ci
        ui <- kronecker(X2, ui)
        ci <- rep(ci, nrow(X2))
        constr <- list(ui = ui, ci = ci)
    } else {
        constr <- do.call(".box_ui_ci", lapply(ret, function(r)
                          attr(r, "constraint")))
    }
    if (length(object) > 1) {
        a <- do.call(".box_char", lapply(ret, function(r) attr(r, "Assign")))
        ret <- do.call(".box", ret)
        attr(ret, "Assign") <- a
    } else {
        ret <- ret[[1]]
    }
    attr(ret, "constraint") <- constr
    return(ret )
}

predict.box_bases <- function(object, newdata, coef, 
                              dim = !is.data.frame(newdata), ...) {

    if (is.logical(dim) & !isTRUE(dim))
        return(predict.basis(object = object, newdata = newdata,
                             coef = coef, dim = dim, ...))

    if (isTRUE(dim)) dim <- sapply(newdata, NROW)

    X <- model.matrix(object = object, data = newdata, 
                      model.matrix = FALSE, dim = dim, ...)
    X <- lapply(X, function(x) as(x, "matrix"))
    X$beta <- array(coef, sapply(X, NCOL))
    lp <- do.call(".cXb", X) 

    vn <- unlist(variable.names(object))
    ### <FIXME> dimensions assessed by first varname; does this work?
    if (length(dim) == 2)
        vn <- vn[vn %in% names(dim)]
    ### </FIXME>
    ret <- .const_array(dim, vn, c(lp))
    ret
}
