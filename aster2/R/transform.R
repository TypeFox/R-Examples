
transformSaturated <- function(parm, data,
    from = c("theta", "phi", "xi", "mu"), to = c("theta", "phi", "xi", "mu"),
    differential, tolerance = 8 * .Machine$double.eps) {
    from <- match.arg(from)
    to <- match.arg(to)
    stopifnot(inherits(data, "asterdata"))
    validasterdata(data)
    stopifnot(is.atomic(parm))
    stopifnot(is.numeric(parm))
    stopifnot(is.finite(parm))
    stopifnot(length(parm) == length(data$repred))
    stopifnot(is.atomic(tolerance))
    stopifnot(is.numeric(tolerance))
    stopifnot(length(tolerance) == 1)
    stopifnot(tolerance > 0)
    fam.set.tolerance(tolerance)
    fam.clear()
    for (i in seq(along = data$families))
        fam.set(data$families[[i]])
    result <- NULL
    if (missing(differential)) {
        if (from == to)
            result <- as.vector(parm)
        if (from == "theta" && to == "phi") {
            if (! is.validThetaNoSetNoClear(data, parm))
                stop("invalid theta vector")
            out <- .C("aster_theta_to_phi",
                nnode = length(parm),
                deriv = as.integer(0),
                pred = as.integer(data$repred),
                group = as.integer(data$regroup),
                code = as.integer(data$recode),
                delta = as.double(data$redelta),
                theta = as.double(parm),
                dtheta = double(0),
                phi = double(length(parm)),
                dphi = double(0),
                PACKAGE = "aster2")
            result <- out$phi
        }
        if (from == "phi" && to == "theta") {
            out <- .C("aster_phi_to_theta",
                nnode = length(parm),
                deriv = as.integer(0),
                pred = as.integer(data$repred),
                group = as.integer(data$regroup),
                code = as.integer(data$recode),
                delta = as.double(data$redelta),
                phi = as.double(parm),
                dphi = double(0),
                theta = double(length(parm)),
                dtheta = double(0),
                PACKAGE = "aster2")
            result <- out$theta
            if (! is.validThetaNoSetNoClear(data, result))
                stop("phi vector maps to invalid theta vector")
        }
        if (from == "theta" && to == "xi") {
            if (! is.validThetaNoSetNoClear(data, parm))
                stop("invalid theta vector")
            out <- .C("aster_theta_to_xi",
                nnode = length(parm),
                deriv = as.integer(0),
                group = as.integer(data$regroup),
                code = as.integer(data$recode),
                delta = as.double(data$redelta),
                theta = as.double(parm),
                dtheta = double(0),
                xi = double(length(parm)),
                dxi = double(0),
                PACKAGE = "aster2")
            result <- out$xi
        }
        if (from == "xi" && to == "theta") {
            if (! is.validXiNoSetNoClear(data, parm))
                stop("invalid xi vector")
            out <- .C("aster_xi_to_theta",
                nnode = length(parm),
                deriv = as.integer(0),
                group = as.integer(data$regroup),
                code = as.integer(data$recode),
                delta = as.double(data$redelta),
                xi = as.double(parm),
                dxi = double(0),
                theta = double(length(parm)),
                dtheta = double(0),
                PACKAGE = "aster2")
            result <- out$theta
        }
        if (from == "xi" && to == "mu") {
            if (! is.validXiNoSetNoClear(data, parm))
                stop("invalid xi vector")
            out <- .C("aster_xi_to_mu",
                nnode = length(parm),
                deriv = as.integer(0),
                pred = as.integer(data$repred),
                initial = as.double(data$initial),
                xi = as.double(parm),
                dxi = double(0),
                mu = double(length(parm)),
                dmu = double(0),
                PACKAGE = "aster2")
            result <- out$mu
        }
        if (from == "mu" && to == "xi") {
            out <- .C("aster_mu_to_xi",
                nnode = length(parm),
                deriv = as.integer(0),
                pred = as.integer(data$repred),
                initial = as.double(data$initial),
                mu = as.double(parm),
                dmu = double(0),
                xi = double(length(parm)),
                dxi = double(0),
                PACKAGE = "aster2")
            if (! is.validXiNoSetNoClear(data, out$xi))
                stop("mu vector maps to invalid xi vector")
            result <- out$xi
        }
        if (from == "theta" && to == "mu") {
            if (! is.validThetaNoSetNoClear(data, parm))
                stop("invalid theta vector")
            out <- .C("aster_theta_to_xi",
                nnode = length(parm),
                deriv = as.integer(0),
                group = as.integer(data$regroup),
                code = as.integer(data$recode),
                delta = as.double(data$redelta),
                theta = as.double(parm),
                dtheta = double(0),
                xi = double(length(parm)),
                dxi = double(0),
                PACKAGE = "aster2")
            out <- .C("aster_xi_to_mu",
                nnode = length(parm),
                deriv = as.integer(0),
                pred = as.integer(data$repred),
                initial = as.double(data$initial),
                xi = out$xi,
                dxi = double(0),
                mu = double(length(parm)),
                dmu = double(0),
                PACKAGE = "aster2")
            result <- out$mu
        }
        if (from == "xi" && to == "phi") {
            if (! is.validXiNoSetNoClear(data, parm))
                stop("invalid xi vector")
            out <- .C("aster_xi_to_theta",
                nnode = length(parm),
                deriv = as.integer(0),
                group = as.integer(data$regroup),
                code = as.integer(data$recode),
                delta = as.double(data$redelta),
                xi = as.double(parm),
                dxi = double(0),
                theta = double(length(parm)),
                dtheta = double(0),
                PACKAGE = "aster2")
            out <- .C("aster_theta_to_phi",
                nnode = length(parm),
                deriv = as.integer(0),
                pred = as.integer(data$repred),
                group = as.integer(data$regroup),
                code = as.integer(data$recode),
                delta = as.double(data$redelta),
                theta = out$theta,
                dtheta = double(0),
                phi = double(length(parm)),
                dphi = double(0),
                PACKAGE = "aster2")
            result <- out$phi
        }
        if (from == "phi" && to == "xi") {
            out <- .C("aster_phi_to_theta",
                nnode = length(parm),
                deriv = as.integer(0),
                pred = as.integer(data$repred),
                group = as.integer(data$regroup),
                code = as.integer(data$recode),
                delta = as.double(data$redelta),
                phi = as.double(parm),
                dphi = double(0),
                theta = double(length(parm)),
                dtheta = double(0),
                PACKAGE = "aster2")
            if (! is.validThetaNoSetNoClear(data, out$theta))
                stop("phi vector maps to invalid theta vector")
            out <- .C("aster_theta_to_xi",
                nnode = length(parm),
                deriv = as.integer(0),
                group = as.integer(data$regroup),
                code = as.integer(data$recode),
                delta = as.double(data$redelta),
                theta = as.double(out$theta),
                dtheta = double(0),
                xi = double(length(parm)),
                dxi = double(0),
                PACKAGE = "aster2")
            result <- out$xi
        }
        if (from == "phi" && to == "mu") {
            out <- .C("aster_phi_to_theta",
                nnode = length(parm),
                deriv = as.integer(0),
                pred = as.integer(data$repred),
                group = as.integer(data$regroup),
                code = as.integer(data$recode),
                delta = as.double(data$redelta),
                phi = as.double(parm),
                dphi = double(0),
                theta = double(length(parm)),
                dtheta = double(0),
                PACKAGE = "aster2")
            if (! is.validThetaNoSetNoClear(data, out$theta))
                stop("phi vector maps to invalid theta vector")
            out <- .C("aster_theta_to_xi",
                nnode = length(parm),
                deriv = as.integer(0),
                group = as.integer(data$regroup),
                code = as.integer(data$recode),
                delta = as.double(data$redelta),
                theta = as.double(out$theta),
                dtheta = double(0),
                xi = double(length(parm)),
                dxi = double(0),
                PACKAGE = "aster2")
            out <- .C("aster_xi_to_mu",
                nnode = length(parm),
                deriv = as.integer(0),
                pred = as.integer(data$repred),
                initial = as.double(data$initial),
                xi = as.double(out$xi),
                dxi = double(0),
                mu = double(length(parm)),
                dmu = double(0),
                PACKAGE = "aster2")
            result <- out$mu
        }
        if (from == "mu" && to == "theta") {
            out <- .C("aster_mu_to_xi",
                nnode = length(parm),
                deriv = as.integer(0),
                pred = as.integer(data$repred),
                initial = as.double(data$initial),
                mu = as.double(parm),
                dmu = double(0),
                xi = double(length(parm)),
                dxi = double(0),
                PACKAGE = "aster2")
            if (! is.validXiNoSetNoClear(data, out$xi))
                stop("mu vector maps to invalid xi vector")
            out <- .C("aster_xi_to_theta",
                nnode = length(parm),
                deriv = as.integer(0),
                group = as.integer(data$regroup),
                code = as.integer(data$recode),
                delta = as.double(data$redelta),
                xi = as.double(out$xi),
                dxi = double(0),
                theta = double(length(parm)),
                dtheta = double(0),
                PACKAGE = "aster2")
            result <- out$theta
        }
        if (from == "mu" && to == "phi") {
            out <- .C("aster_mu_to_xi",
                nnode = length(parm),
                deriv = as.integer(0),
                pred = as.integer(data$repred),
                initial = as.double(data$initial),
                mu = as.double(parm),
                dmu = double(0),
                xi = double(length(parm)),
                dxi = double(0),
                PACKAGE = "aster2")
            if (! is.validXiNoSetNoClear(data, out$xi))
                stop("mu vector maps to invalid xi vector")
            out <- .C("aster_xi_to_theta",
                nnode = length(parm),
                deriv = as.integer(0),
                group = as.integer(data$regroup),
                code = as.integer(data$recode),
                delta = as.double(data$redelta),
                xi = as.double(out$xi),
                dxi = double(0),
                theta = double(length(parm)),
                dtheta = double(0),
                PACKAGE = "aster2")
            out <- .C("aster_theta_to_phi",
                nnode = length(parm),
                deriv = as.integer(0),
                pred = as.integer(data$repred),
                group = as.integer(data$regroup),
                code = as.integer(data$recode),
                delta = as.double(data$redelta),
                theta = as.double(out$theta),
                dtheta = double(0),
                phi = double(length(parm)),
                dphi = double(0),
                PACKAGE = "aster2")
            result <- out$phi
        }
    } else {
        ##### ! missing(differential) #####
        stopifnot(is.atomic(differential))
        stopifnot(is.numeric(differential))
        stopifnot(is.finite(differential))
        stopifnot(length(differential) == length(parm))
        if (from == to)
            result <- as.vector(differential)
        if (from == "theta" && to == "phi") {
            if (! is.validThetaNoSetNoClear(data, parm))
                stop("invalid theta vector")
            out <- .C("aster_theta_to_phi",
                nnode = length(parm),
                deriv = as.integer(1),
                pred = as.integer(data$repred),
                group = as.integer(data$regroup),
                code = as.integer(data$recode),
                delta = as.double(data$redelta),
                theta = as.double(parm),
                dtheta = as.double(differential),
                phi = double(length(parm)),
                dphi = double(length(parm)),
                PACKAGE = "aster2")
            result <- out$dphi
        }
        if (from == "phi" && to == "theta") {
            out <- .C("aster_phi_to_theta",
                nnode = length(parm),
                deriv = as.integer(1),
                pred = as.integer(data$repred),
                group = as.integer(data$regroup),
                code = as.integer(data$recode),
                delta = as.double(data$redelta),
                phi = as.double(parm),
                dphi = as.double(differential),
                theta = double(length(parm)),
                dtheta = double(length(parm)),
                PACKAGE = "aster2")
            result <- out$dtheta
            if (! is.validThetaNoSetNoClear(data, out$theta))
                stop("phi vector maps to invalid theta vector")
        }
        if (from == "theta" && to == "xi") {
            if (! is.validThetaNoSetNoClear(data, parm))
                stop("invalid theta vector")
            out <- .C("aster_theta_to_xi",
                nnode = length(parm),
                deriv = as.integer(1),
                group = as.integer(data$regroup),
                code = as.integer(data$recode),
                delta = as.double(data$redelta),
                theta = as.double(parm),
                dtheta = as.double(differential),
                xi = double(length(parm)),
                dxi = double(length(parm)),
                PACKAGE = "aster2")
            result <- out$dxi
        }
        if (from == "xi" && to == "theta") {
            if (! is.validXiNoSetNoClear(data, parm))
                stop("invalid xi vector")
            out <- .C("aster_xi_to_theta",
                nnode = length(parm),
                deriv = as.integer(1),
                group = as.integer(data$regroup),
                code = as.integer(data$recode),
                delta = as.double(data$redelta),
                xi = as.double(parm),
                dxi = as.double(differential),
                theta = double(length(parm)),
                dtheta = double(length(parm)),
                PACKAGE = "aster2")
            result <- out$dtheta
        }
        if (from == "xi" && to == "mu") {
            if (! is.validXiNoSetNoClear(data, parm))
                stop("invalid xi vector")
            out <- .C("aster_xi_to_mu",
                nnode = length(parm),
                deriv = as.integer(1),
                pred = as.integer(data$repred),
                initial = as.double(data$initial),
                xi = as.double(parm),
                dxi = as.double(differential),
                mu = double(length(parm)),
                dmu = double(length(parm)),
                PACKAGE = "aster2")
            result <- out$dmu
        }
        if (from == "mu" && to == "xi") {
            out <- .C("aster_mu_to_xi",
                nnode = length(parm),
                deriv = as.integer(1),
                pred = as.integer(data$repred),
                initial = as.double(data$initial),
                mu = as.double(parm),
                dmu = as.double(differential),
                xi = double(length(parm)),
                dxi = double(length(parm)),
                PACKAGE = "aster2")
            if (! is.validXiNoSetNoClear(data, out$xi))
                stop("mu vector maps to invalid xi vector")
            result <- out$dxi
        }
        if (from == "theta" && to == "mu") {
            if (! is.validThetaNoSetNoClear(data, parm))
                stop("invalid theta vector")
            out <- .C("aster_theta_to_xi",
                nnode = length(parm),
                deriv = as.integer(1),
                group = as.integer(data$regroup),
                code = as.integer(data$recode),
                delta = as.double(data$redelta),
                theta = as.double(parm),
                dtheta = as.double(differential),
                xi = double(length(parm)),
                dxi = double(length(parm)),
                PACKAGE = "aster2")
            out <- .C("aster_xi_to_mu",
                nnode = length(parm),
                deriv = as.integer(1),
                pred = as.integer(data$repred),
                initial = as.double(data$initial),
                xi = out$xi,
                dxi = out$dxi,
                mu = double(length(parm)),
                dmu = double(length(parm)),
                PACKAGE = "aster2")
            result <- out$dmu
        }
        if (from == "xi" && to == "phi") {
            if (! is.validXiNoSetNoClear(data, parm))
                stop("invalid xi vector")
            out <- .C("aster_xi_to_theta",
                nnode = length(parm),
                deriv = as.integer(1),
                group = as.integer(data$regroup),
                code = as.integer(data$recode),
                delta = as.double(data$redelta),
                xi = as.double(parm),
                dxi = as.double(differential),
                theta = double(length(parm)),
                dtheta = double(length(parm)),
                PACKAGE = "aster2")
            out <- .C("aster_theta_to_phi",
                nnode = length(parm),
                deriv = as.integer(1),
                pred = as.integer(data$repred),
                group = as.integer(data$regroup),
                code = as.integer(data$recode),
                delta = as.double(data$redelta),
                theta = out$theta,
                dtheta = out$dtheta,
                phi = double(length(parm)),
                dphi = double(length(parm)),
                PACKAGE = "aster2")
            result <- out$dphi
        }
        if (from == "phi" && to == "xi") {
            out <- .C("aster_phi_to_theta",
                nnode = length(parm),
                deriv = as.integer(1),
                pred = as.integer(data$repred),
                group = as.integer(data$regroup),
                code = as.integer(data$recode),
                delta = as.double(data$redelta),
                phi = as.double(parm),
                dphi = as.double(differential),
                theta = double(length(parm)),
                dtheta = double(length(parm)),
                PACKAGE = "aster2")
            if (! is.validThetaNoSetNoClear(data, out$theta))
                stop("phi vector maps to invalid theta vector")
            out <- .C("aster_theta_to_xi",
                nnode = length(parm),
                deriv = as.integer(1),
                group = as.integer(data$regroup),
                code = as.integer(data$recode),
                delta = as.double(data$redelta),
                theta = as.double(out$theta),
                dtheta = as.double(out$dtheta),
                xi = double(length(parm)),
                dxi = double(length(parm)),
                PACKAGE = "aster2")
            result <- out$dxi
        }
        if (from == "phi" && to == "mu") {
            out <- .C("aster_phi_to_theta",
                nnode = length(parm),
                deriv = as.integer(1),
                pred = as.integer(data$repred),
                group = as.integer(data$regroup),
                code = as.integer(data$recode),
                delta = as.double(data$redelta),
                phi = as.double(parm),
                dphi = as.double(differential),
                theta = double(length(parm)),
                dtheta = double(length(parm)),
                PACKAGE = "aster2")
            if (! is.validThetaNoSetNoClear(data, out$theta))
                stop("phi vector maps to invalid theta vector")
            out <- .C("aster_theta_to_xi",
                nnode = length(parm),
                deriv = as.integer(1),
                group = as.integer(data$regroup),
                code = as.integer(data$recode),
                delta = as.double(data$redelta),
                theta = as.double(out$theta),
                dtheta = as.double(out$dtheta),
                xi = double(length(parm)),
                dxi = double(length(parm)),
                PACKAGE = "aster2")
            out <- .C("aster_xi_to_mu",
                nnode = length(parm),
                deriv = as.integer(1),
                pred = as.integer(data$repred),
                initial = as.double(data$initial),
                xi = as.double(out$xi),
                dxi = as.double(out$dxi),
                mu = double(length(parm)),
                dmu = double(length(parm)),
                PACKAGE = "aster2")
            result <- out$dmu
        }
        if (from == "mu" && to == "theta") {
            out <- .C("aster_mu_to_xi",
                nnode = length(parm),
                deriv = as.integer(1),
                pred = as.integer(data$repred),
                initial = as.double(data$initial),
                mu = as.double(parm),
                dmu = as.double(differential),
                xi = double(length(parm)),
                dxi = double(length(parm)),
                PACKAGE = "aster2")
            if (! is.validXiNoSetNoClear(data, out$xi))
                stop("mu vector maps to invalid xi vector")
            out <- .C("aster_xi_to_theta",
                nnode = length(parm),
                deriv = as.integer(1),
                group = as.integer(data$regroup),
                code = as.integer(data$recode),
                delta = as.double(data$redelta),
                xi = as.double(out$xi),
                dxi = as.double(out$dxi),
                theta = double(length(parm)),
                dtheta = double(length(parm)),
                PACKAGE = "aster2")
            result <- out$dtheta
        }
        if (from == "mu" && to == "phi") {
            out <- .C("aster_mu_to_xi",
                nnode = length(parm),
                deriv = as.integer(1),
                pred = as.integer(data$repred),
                initial = as.double(data$initial),
                mu = as.double(parm),
                dmu = as.double(differential),
                xi = double(length(parm)),
                dxi = double(length(parm)),
                PACKAGE = "aster2")
            if (! is.validXiNoSetNoClear(data, out$xi))
                stop("mu vector maps to invalid xi vector")
            out <- .C("aster_xi_to_theta",
                nnode = length(parm),
                deriv = as.integer(1),
                group = as.integer(data$regroup),
                code = as.integer(data$recode),
                delta = as.double(data$redelta),
                xi = as.double(out$xi),
                dxi = as.double(out$dxi),
                theta = double(length(parm)),
                dtheta = double(length(parm)),
                PACKAGE = "aster2")
            out <- .C("aster_theta_to_phi",
                nnode = length(parm),
                deriv = as.integer(1),
                pred = as.integer(data$repred),
                group = as.integer(data$regroup),
                code = as.integer(data$recode),
                delta = as.double(data$redelta),
                theta = as.double(out$theta),
                dtheta = as.double(out$dtheta),
                phi = double(length(parm)),
                dphi = double(length(parm)),
                PACKAGE = "aster2")
            result <- out$dphi
        }
    }
    fam.reset.tolerance()
    fam.clear()
    if (is.null(result)) stop("Can't happen, this is a bug")
    return(result)
}

transformConditional <- function(parm, modmat, data, from = "beta",
    to = c("theta", "phi", "xi", "mu"), differential, offset,
    tolerance = 8 * .Machine$double.eps) {
    from <- match.arg(from)
    to <- match.arg(to)
    stopifnot(inherits(data, "asterdata"))
    validasterdata(data)
    stopifnot(is.atomic(parm))
    stopifnot(is.numeric(parm))
    stopifnot(is.finite(parm))
    stopifnot(is.atomic(modmat))
    stopifnot(is.numeric(modmat))
    stopifnot(is.finite(modmat))
    stopifnot(is.matrix(modmat))
    stopifnot(nrow(modmat) == length(data$repred))
    stopifnot(ncol(modmat) == length(parm))
    stopifnot(is.atomic(tolerance))
    stopifnot(is.numeric(tolerance))
    stopifnot(length(tolerance) == 1)
    stopifnot(tolerance > 0)
    if (missing(offset))
        offset <- rep(0, length(data$repred))
    stopifnot(is.atomic(offset))
    stopifnot(is.numeric(offset))
    stopifnot(is.finite(offset))
    stopifnot(length(offset) == length(data$repred))
    result <- NULL
    if (missing(differential)) {
        result <- transformSaturated(offset + modmat %*% parm, data,
            from = "theta", to = to, tolerance = tolerance)
    } else {
        stopifnot(is.atomic(differential))
        stopifnot(is.numeric(differential))
        stopifnot(is.finite(differential))
        stopifnot(length(differential) == length(parm))
        result <- transformSaturated(offset + modmat %*% parm, data,
            from = "theta", to = to, differential = modmat %*% differential,
            tolerance = tolerance)
    }
    return(result)
}

transformUnconditional <- function(parm, modmat, data,
    from = c("beta", "tau"), to = c("beta", "theta", "phi", "xi", "mu", "tau"),
    differential, offset, tolerance = 8 * .Machine$double.eps) {
    from <- match.arg(from)
    to <- match.arg(to)
    stopifnot(inherits(data, "asterdata"))
    validasterdata(data)
    stopifnot(is.atomic(parm))
    stopifnot(is.numeric(parm))
    stopifnot(is.finite(parm))
    stopifnot(is.atomic(modmat))
    stopifnot(is.numeric(modmat))
    stopifnot(is.finite(modmat))
    stopifnot(is.matrix(modmat))
    stopifnot(nrow(modmat) == length(data$repred))
    stopifnot(ncol(modmat) == length(parm))
    stopifnot(is.atomic(tolerance))
    stopifnot(is.numeric(tolerance))
    stopifnot(length(tolerance) == 1)
    stopifnot(tolerance > 0)
    if (missing(offset))
        offset <- rep(0, length(data$repred))
    stopifnot(is.atomic(offset))
    stopifnot(is.numeric(offset))
    stopifnot(is.finite(offset))
    stopifnot(length(offset) == length(data$repred))
    result <- NULL
    if (missing(differential)) {
        if (from == to)
            result <- parm
        if (from == "beta" && to %in% c("theta", "phi", "xi", "mu"))
            result <- transformSaturated(offset + modmat %*% parm, data,
                from = "phi", to = to, tolerance = tolerance)
        if (from == "beta" && to == "tau")
            result <- as.vector(t(modmat) %*%
                transformSaturated(offset + modmat %*% parm, data,
                from = "phi", to = "mu", tolerance = tolerance))
    } else {
        stopifnot(is.atomic(differential))
        stopifnot(is.numeric(differential))
        stopifnot(is.finite(differential))
        stopifnot(length(differential) == length(parm))
        if (from == to)
            result <- differential
        if (from == "beta" && to %in% c("theta", "phi", "xi", "mu"))
            result <- transformSaturated(offset + modmat %*% parm, data,
                from = "phi", to = to, differential = modmat %*% differential,
                tolerance = tolerance)
        if (from == "beta" && to == "tau")
            result <- as.vector(t(modmat) %*%
                transformSaturated(offset + modmat %*% parm, data,
                from = "phi", to = "mu",
                differential = modmat %*% differential, tolerance = tolerance))
    }
    if (from == "tau" & to != "tau") {
        consmat <- constancy(data, parm.type = "phi")
        theta.start <- starting(data)
        phi.start <- transformSaturated(theta.start, data,
            from = "theta", to = "phi")
        modmat.augmented <- cbind(t(consmat), modmat)
        foo <- qr(modmat.augmented)
        bar <- qr.coef(foo, phi.start - offset)
        bar <- bar[seq(nrow(consmat) + 1, length(bar))]
        baz <- (! is.na(bar))
        beta.start <- bar
        beta.start[is.na(beta.start)] <- 0
        phi.start <- as.numeric(offset + modmat %*% beta.start)
        theta.start <- transformSaturated(phi.start, data,
            from = "phi", to = "theta", tolerance = tolerance)
        if (! is.validtheta(data, theta.start)) {
            stop(paste("cannot find valid beta to start iteration,",
               " try different offset"))
        }
        my.modmat <- modmat[ , baz]
        my.tau <- parm[baz]
        my.beta.start <- beta.start[baz]
        gradfun <- function(x) my.tau - transformUnconditional(x,
            my.modmat, data, from = "beta", to = "tau",
            offset = offset)
        hessfun <- function(x) - jacobian(x, data,
            transform = "unconditional",
            from = "beta", to = "tau", modmat = my.modmat,
            offset = offset)
        out <- outer.loop(my.beta.start, gradfun, hessfun, my.tau,
            fudge = 0.05, tol1 = sqrt(.Machine$double.eps), tol2 = 0.10)
        my.beta <- out$x
        beta <- rep(0, length(parm))
        beta[baz] <- my.beta
        if (missing(differential)) {
            if (to == "beta") {
                result <- beta
            } else {
                phi <- as.numeric(offset + modmat %*% beta)
                if (to == "phi") {
                    result <- phi
                } else {
                    result <- transformSaturated(phi, data,
                        from = "phi", to = to, tolerance = tolerance)
                }
            }
        } else {
            my.dtau <- differential[baz]
            my.hess <- hessfun(my.beta)
            my.dbeta <- solve(- my.hess, my.dtau)
            dbeta <- rep(0, length(parm))
            dbeta[baz] <- my.dbeta
            if (to == "beta") {
                result <- dbeta
            } else {
                phi <- as.numeric(offset + modmat %*% beta)
                dphi <- as.numeric(modmat %*% dbeta)
                if (to == "phi") {
                    result <- dphi
                } else {
                    result <- transformSaturated(phi, data,
                        from = "phi", to = to, differential = dphi,
                        tolerance = tolerance)
                }
            }
        }
    }
    if (is.null(result)) stop("transformation not implemented yet")
    return(result)
}

jacobian <- function(parm, data,
    transform = c("saturated", "conditional", "unconditional"),
    from = c("beta", "theta", "phi", "xi", "mu", "tau"),
    to = c("beta", "theta", "phi", "xi", "mu", "tau"),
    modmat, offset, tolerance = 8 * .Machine$double.eps) {
    transform <- match.arg(transform)
    from <- match.arg(from)
    to <- match.arg(to)
    foo <- switch(transform, saturated = transformSaturated,
        conditional = transformConditional,
        unconditional = transformUnconditional)
    fred <- as.list(args(foo))
    fred.from <- eval(fred$from)
    fred.to <- eval(fred$to)
    if (! (from %in% fred.from))
        stop(paste("from = \"", from, "\" not valid for transform = \"",
            transform, "\"", sep = ""))
    if (! (to %in% fred.to))
        stop(paste("to = \"", to, "\" not valid for transform = \"",
            transform, "\"", sep = ""))
    stopifnot(inherits(data, "asterdata"))
    validasterdata(data)
    stopifnot(is.atomic(parm))
    stopifnot(is.numeric(parm))
    stopifnot(is.finite(parm))
    stopifnot(is.atomic(tolerance))
    stopifnot(is.numeric(tolerance))
    stopifnot(length(tolerance) == 1)
    stopifnot(tolerance > 0)
    if (transform == "saturated") {
        stopifnot(length(parm) == length(data$repred))
        bar <- function(baz) foo(parm, data, from = from, to = to,
            differential = baz, tolerance = tolerance)
    } else {
        stopifnot(is.atomic(modmat))
        stopifnot(is.numeric(modmat))
        stopifnot(is.finite(modmat))
        stopifnot(is.matrix(modmat))
        stopifnot(nrow(modmat) == length(data$repred))
        stopifnot(ncol(modmat) == length(parm))
        if (missing(offset))
            offset <- rep(0, length(data$repred))
        stopifnot(is.atomic(offset))
        stopifnot(is.numeric(offset))
        stopifnot(is.finite(offset))
        stopifnot(length(offset) == length(data$repred))
        bar <- function(baz) foo(parm, modmat, data, from = from, to = to,
            differential = baz, offset = offset, tolerance = tolerance)
    }
    zeros <- rep(0, length(parm))
    e1 <- zeros
    e1[1] <- 1
    v1 <- bar(e1)
    result <- matrix(NA, nrow = length(v1), ncol = length(e1))
    result[ , 1] <- v1
    if (length(parm) > 1) {
        for (j in 2:length(parm)) {
            ej <- zeros
            ej[j] <- 1
            vj <- bar(ej)
            result[ , j] <- vj
        }
    }
    return(result)
}

validtheta <- function(data, theta, tolerance = 8 * .Machine$double.eps) {
    stopifnot(inherits(data, "asterdata"))
    validasterdata(data)
    stopifnot(is.atomic(theta))
    stopifnot(is.numeric(theta))
    stopifnot(is.finite(theta))
    stopifnot(length(theta) == length(data$repred))
    stopifnot(is.atomic(tolerance))
    stopifnot(is.numeric(tolerance))
    stopifnot(length(tolerance) == 1)
    stopifnot(tolerance > 0)
    fam.set.tolerance(tolerance)
    fam.clear()
    for (i in seq(along = data$families))
        fam.set(data$families[[i]])
    out <- .C("aster_validate_theta", nnode = length(theta),
        group = as.integer(data$regroup), code = as.integer(data$recode),
        delta = as.double(data$redelta), theta = as.double(theta),
        PACKAGE = "aster2")
    fam.reset.tolerance()
    fam.clear()
    invisible(TRUE)
}

is.validtheta <- function(data, theta, tolerance = 8 * .Machine$double.eps) {
    out <- try(validtheta(data, theta, tolerance), silent = TRUE)
    return(! inherits(out, "try-error"))
}

## also need for internal use only an is.validtheta that doesn't set families
## and doesn't clear them either

is.validThetaNoSetNoClear <- function(data, theta) {
    stopifnot(inherits(data, "asterdata"))
    stopifnot(is.atomic(theta))
    stopifnot(is.numeric(theta))
    stopifnot(is.finite(theta))
    stopifnot(length(theta) == length(data$repred))
    out <- try(.C("aster_validate_theta", nnode = length(theta),
        group = as.integer(data$regroup), code = as.integer(data$recode),
        delta = as.double(data$redelta), theta = as.double(theta),
        PACKAGE = "aster2"), silent = TRUE)
    return(! inherits(out, "try-error"))
}

validxi <- function(data, xi, tolerance = 8 * .Machine$double.eps) {
    stopifnot(inherits(data, "asterdata"))
    validasterdata(data)
    stopifnot(is.atomic(xi))
    stopifnot(is.numeric(xi))
    stopifnot(is.finite(xi))
    stopifnot(length(xi) == length(data$repred))
    stopifnot(is.atomic(tolerance))
    stopifnot(is.numeric(tolerance))
    stopifnot(length(tolerance) == 1)
    stopifnot(tolerance > 0)
    fam.set.tolerance(tolerance)
    fam.clear()
    for (i in seq(along = data$families))
        fam.set(data$families[[i]])
    out <- .C("aster_validate_xi", nnode = length(xi),
        group = as.integer(data$regroup), code = as.integer(data$recode),
        delta = as.double(data$redelta), xi = as.double(xi),
        PACKAGE = "aster2")
    fam.reset.tolerance()
    fam.clear()
    invisible(TRUE)
}

is.validxi <- function(data, xi, tolerance = 8 * .Machine$double.eps) {
    stopifnot(inherits(data, "asterdata"))
    out <- try(validxi(data, xi, tolerance), silent = TRUE)
    return(! inherits(out, "try-error"))
}

## also need for internal use only an is.validtheta that doesn't set families
## and doesn't clear them either

is.validXiNoSetNoClear <- function(data, xi) {
    stopifnot(inherits(data, "asterdata"))
    stopifnot(is.atomic(xi))
    stopifnot(is.numeric(xi))
    stopifnot(is.finite(xi))
    stopifnot(length(xi) == length(data$repred))
    out <- try(.C("aster_validate_xi", nnode = length(xi),
        group = as.integer(data$regroup), code = as.integer(data$recode),
        delta = as.double(data$redelta), xi = as.double(xi),
        PACKAGE = "aster2"), silent = TRUE)
    return(! inherits(out, "try-error"))
}

