twostage <- function(time, delta, group, nboot, alpha = 0.05, eps = 0.1) {

    stopifnot(is.numeric(time))
    stopifnot(is.numeric(delta))
    stopifnot(is.numeric(group))
    stopifnot(is.numeric(nboot))
    stopifnot(is.numeric(alpha))
    stopifnot(is.numeric(eps))

    n <- length(time)
    stopifnot(length(delta) == length(time))
    stopifnot(length(group) == length(time))
    stopifnot(length(nboot) == 1)
    stopifnot(length(alpha) == 1)
    stopifnot(length(eps) == 1)

    stopifnot(all(time >= 0))
    stopifnot(all(delta %in% c(0, 1)))
    stopifnot(all(group %in% c(0, 1)))
    stopifnot(nboot == as.integer(nboot))
    stopifnot(nboot > 0)
    stopifnot(0 < alpha & alpha < 1)
    stopifnot(0 < eps & eps < 1)

    out <- .Fortran("TWOSTAGE",
        WORK = as.integer(n),
        DATASIZE = as.integer(n),
        T = as.double(time),
        DELTA = as.integer(delta),
        GP = as.integer(group),
        BOOTSN = as.integer(nboot),
        ALPHA = as.double(alpha),
        EPS = as.double(eps),
        LRPV = double(1),
        MTPV = double(1),
        TSPV = double(1),
        PACKAGE = "TSHRC")

    result <- c(out$LRPV, out$MTPV, out$TSPV)
    names(result) <- c("LRPV", "MTPV", "TSPV")
    return(result)
}

