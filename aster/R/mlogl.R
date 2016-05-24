
mlogl <- function(parm, pred, fam, x, root, modmat, deriv = 0,
    type = c("unconditional", "conditional"), famlist = fam.default(),
    origin, origin.type = c("model.type", "unconditional", "conditional"))
{
    type <- match.arg(type)
    origin.type <- match.arg(origin.type)

    stopifnot(is.element(fam, seq(along = famlist)))
    setfam(famlist)

    ##### should probably be separate function not snarf-and-barf from aster
    nnode <- length(fam)
    stopifnot(is.numeric(x))
    nind <- floor(length(x) / nnode)
    stopifnot(length(x) == nind * nnode)
    if (missing(origin)) {
        origin <- .C("aster_default_origin",
            nind = as.integer(nind),
            nnode = as.integer(nnode),
            fam = as.integer(fam),
            theta = matrix(as.double(0), nind, nnode),
            PACKAGE = "aster")$theta
        if (type == "unconditional")
            origin <- .C("aster_theta2phi",
                nind = as.integer(nind),
                nnode = as.integer(nnode),
                pred = as.integer(pred),
                fam = as.integer(fam),
                theta = origin,
                phi = matrix(as.double(0), nind, nnode),
                PACKAGE = "aster")$phi
    } else {
        stopifnot(is.numeric(origin))
        storage.mode(origin) <- "double"
        stopifnot(length(origin) == nind * nnode)
        stopifnot(all(is.finite(origin)))
        if (is.matrix(origin)) {
            stopifnot(nrow(origin) == nind)
            stopifnot(ncol(origin) == nnode)
        } else {
            origin <- matrix(origin, nind, nnode)
        }
        if (origin.type == "model.type")
            origin.type <- type
        if (type == "unconditional" && origin.type == "conditional")
            origin <- .C("aster_theta2phi",
                nind = as.integer(nind),
                nnode = as.integer(nnode),
                pred = as.integer(pred),
                fam = as.integer(fam),
                theta = origin,
                phi = matrix(as.double(0), nind, nnode),
                PACKAGE = "aster")$phi
        if (type == "conditional" && origin.type == "unconditional")
            origin <- .C("aster_phi2theta",
                nind = as.integer(nind),
                nnode = as.integer(nnode),
                pred = as.integer(pred),
                fam = as.integer(fam),
                phi = origin,
                theta = matrix(as.double(0), nind, nnode),
                PACKAGE = "aster")$theta
    }

    result <- mloglhelper(parm, pred, fam, x, root, modmat, origin,
        deriv, type)
    clearfam()
    return(result)
}

mloglhelper <- function(parm, pred, fam, x, root, modmat, origin,
    deriv = 0, type = c("unconditional", "conditional"))
{
    type <- match.arg(type)

    stopifnot(is.numeric(x))
    stopifnot(is.numeric(root))
    stopifnot(is.numeric(modmat))
    stopifnot(is.numeric(parm))
    stopifnot(is.numeric(pred))
    stopifnot(is.numeric(fam))
    stopifnot(is.numeric(deriv))

    stopifnot(all(pred == as.integer(pred)))
    stopifnot(all(fam == as.integer(fam)))
    stopifnot(all(deriv == as.integer(deriv)))

    stopifnot(all(pred < seq(along = pred)))
    stopifnot(is.element(deriv, 0:2))
    stopifnot(length(deriv) == 1)

    ncoef <- length(parm)
    nnode <- length(pred)

    foo <- dim(modmat)
    if (ncoef != foo[length(foo)])
        stop("last dimension of modmat not length(parm)")
    if (length(x) != prod(foo[- length(foo)]))
        stop("product of all but last dimension of modmat not length(x)")
    stopifnot(length(x) == length(root))
    stopifnot(length(pred) == length(fam))

    nind <- length(x) / nnode
    if(nind != as.integer(nind))
        stop("length(x) not multiple of length(pred)")

    stopifnot(length(origin) == nind * nnode)

    if (type == "unconditional") {
        out <- .C("aster_mlogl_unco",
            nind = as.integer(nind),
            nnode = as.integer(nnode),
            ncoef = as.integer(ncoef),
            pred = as.integer(pred),
            fam = as.integer(fam),
            deriv = as.integer(deriv),
            beta = as.double(parm),
            root = as.double(root),
            x = as.double(x),
            origin = as.double(origin),
            modmat = as.double(modmat),
            value = double(1),
            gradient = double(ncoef),
            hessian = matrix(as.double(0), ncoef, ncoef),
            PACKAGE = "aster")
    } else {
        out <- .C("aster_mlogl_cond",
            nind = as.integer(nind),
            nnode = as.integer(nnode),
            ncoef = as.integer(ncoef),
            pred = as.integer(pred),
            fam = as.integer(fam),
            deriv = as.integer(deriv),
            beta = as.double(parm),
            root = as.double(root),
            x = as.double(x),
            origin = as.double(origin),
            modmat = as.double(modmat),
            value = double(1),
            gradient = double(ncoef),
            hessian = matrix(as.double(0), ncoef, ncoef),
            PACKAGE = "aster")
    }

    if (out$value == Inf) {
        out$gradient <- NaN * out$gradient
        out$hessian <- NaN * out$hessian
    }

    if (deriv == 0)
        return(list(value = out$value))
    if (deriv == 1)
        return(list(value = out$value, gradient = out$gradient))
    if (deriv == 2)
        return(list(value = out$value, gradient = out$gradient,
            hessian = out$hessian))
}
