
constancy <- function(data, parm.type = c("theta", "phi")) {
    parm.type <- match.arg(parm.type)
    stopifnot(inherits(data, "asterdata"))
    validasterdata(data)
    fam.clear()
    for (i in seq(along = data$families))
        fam.set(data$families[[i]])
    result <- .Call("aster_constancy",
        as.integer(data$repred),
        as.integer(data$regroup),
        as.integer(data$recode),
        as.double(data$redelta),
        parm.type == "theta",
        PACKAGE = "aster2")
    fam.clear()
    return(result)
}

is.same <- function(parm1, parm2, data, parm.type = c("theta", "phi"),
    tolerance = sqrt(.Machine$double.eps)) {
    parm.type <- match.arg(parm.type)
    stopifnot(is.atomic(parm1))
    stopifnot(is.numeric(parm1))
    stopifnot(all(is.finite(parm1)))
    stopifnot(is.atomic(parm2))
    stopifnot(is.numeric(parm2))
    stopifnot(all(is.finite(parm2)))
    stopifnot(inherits(data, "asterdata"))
    validasterdata(data)
    stopifnot(length(parm1) == length(data$repred))
    stopifnot(length(parm2) == length(data$repred))
    stopifnot(is.atomic(tolerance))
    stopifnot(is.numeric(tolerance))
    stopifnot(length(tolerance) == 1)
    stopifnot(tolerance > 0)

    cmat <- constancy(data, parm.type = parm.type)
    foo <- qr(t(cmat), lapack = TRUE)
    bar <- qr.resid(foo, parm1 - parm2)
    return(all(abs(bar) < tolerance))
}

