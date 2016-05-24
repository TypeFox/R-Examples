
# a1 and b1 define inequality constraints (a1 %*% x <= b1)
# a2 and b2 define inequality constraints (a2 %*% x == b2)
# as in the makeH function in the rcdd package
# otherwise arguments are like metrop in the mcmc package
# we add to a1 and b1 the unit simplex constraints (all(x >= 0))
# we add to a2 and b2 the unit simplex constraints (sum(x) == 1)
# a1 and b1 can be missing
# a2 and b2 can be missing

# this file is now revised to follow the design document hitrun3.Rnw
# (in the devel directory) rather than the older hitrun.Rnw and hitrun2.Rnw

hitrun <- function(alpha, ...)
    UseMethod("hitrun")

hitrun.hitrun <- function(alpha, nbatch, blen, nspac, outmat, debug, ...)
{
    if (missing(nbatch)) nbatch <- alpha$nbatch
    if (missing(blen)) blen <- alpha$blen
    if (missing(nspac)) nspac <- alpha$nspac
    if (missing(outmat)) outmat <- alpha$outmat
    if (missing(debug)) debug <- alpha$debug

    assign(".Random.seed", alpha$final.seed, .GlobalEnv)
    basis <- alpha$basis
    origin <- alpha$origin
    amat <- alpha$amat
    bvec <- alpha$bvec

    out <- hitrunHelper(alpha$alpha, alpha$final, nbatch, blen, nspac,
        origin, basis, amat, bvec, outmat, debug)

    class(out) <- "hitrun"
    return(out)
}

hitrun.default <- function(alpha, a1 = NULL, b1 = NULL, a2 = NULL, b2 = NULL, 
     nbatch = 1, blen = 1, nspac = 1, outmat = NULL, debug = FALSE,
     stop.if.implied.equalities = FALSE, ...)
{
    if (! exists(".Random.seed")) runif(1)
    saveseed <- .Random.seed

    stopifnot(is.numeric(alpha))
    stopifnot(is.finite(alpha))
    stopifnot(is.vector(alpha))
    stopifnot(length(alpha) >= 2)
    stopifnot(alpha > 0)
    stopifnot(is.null(a1) == is.null(b1))
    if (! is.null(a1)) {
        stopifnot(is.numeric(a1) || is.character(a1))
        stopifnot(is.matrix(a1))
        if (is.numeric(a1)) {
            stopifnot(is.finite(a1))
            a1 <- d2q(a1)
        } else {
            a1 <- try(q2q(a1))
            if (inherits(a1, "try-error"))
            stop("'a1' character but not GMP rational")
        }
        stopifnot(is.numeric(b1) || is.character(b1))
        stopifnot(is.vector(b1))
        if (is.numeric(b1)) {
            stopifnot(is.finite(b1))
            b1 <- d2q(b1)
        } else {
            b1 <- try(q2q(b1))
            if (inherits(b1, "try-error"))
            stop("'b1' character but not GMP rational")
        }
        stopifnot(nrow(a1) == length(b1))
        stopifnot(ncol(a1) == length(alpha))
    }
    stopifnot(is.null(a2) == is.null(b2))
    if (! is.null(a2)) {
        stopifnot(is.numeric(a2) || is.character(a2))
        stopifnot(is.matrix(a2))
        if (is.numeric(a2)) {
            stopifnot(is.finite(a2))
            a2 <- d2q(a2)
        } else {
            a2 <- try(q2q(a2))
            if (inherits(a2, "try-error"))
            stop("'a2' character but not GMP rational")
        }
        stopifnot(is.numeric(b2) || is.character(b2))
        stopifnot(is.vector(b2))
        if (is.numeric(b2)) {
            stopifnot(is.finite(b2))
            b2 <- d2q(b2)
        } else {
            b2 <- try(q2q(b2))
            if (inherits(b2, "try-error"))
            stop("'b2' character but not GMP rational")
        }
        stopifnot(nrow(a2) == length(b2))
        stopifnot(ncol(a2) == length(alpha))
    }

    if (! is.null(outmat)) {
        stopifnot(is.numeric(outmat))
        stopifnot(is.finite(outmat))
        stopifnot(is.matrix(outmat))
        stopifnot(ncol(outmat) == length(alpha))
    }

    stopifnot(is.numeric(nbatch))
    stopifnot(nbatch == as.integer(nbatch))
    stopifnot(length(nbatch) == 1)
    stopifnot(nbatch >= 1)
    stopifnot(is.numeric(blen))
    stopifnot(blen == as.integer(blen))
    stopifnot(length(blen) == 1)
    stopifnot(blen >= 1)
    stopifnot(is.numeric(nspac))
    stopifnot(nspac == as.integer(nspac))
    stopifnot(length(nspac) == 1)
    stopifnot(nspac >= 1)
    stopifnot(is.logical(debug))
    stopifnot(length(debug) == 1)
    stopifnot(is.logical(stop.if.implied.equalities))
    stopifnot(length(stop.if.implied.equalities) == 1)

    d <- length(alpha)
    # add unit simplex constraints
    a1 <- rbind(a1, d2q(- diag(d)))
    b1 <- c(b1, rep("0", d))
    a2 <- rbind(a2, rep("1", d))
    b2 <- c(b2, "1")

    hrep1 <- makeH(a1, b1, a2, b2)

    # step 1: discover whether linearity step is necessary.

    start.time.lin <- proc.time()

    hrep6 <- cbind(hrep1, "0")
    hrep6[hrep6[ , 1] == "0", length(alpha) + 3] <- "-1"
    grad6 <- c(rep("0", d), "1")
    lout <- lpcdd(hrep6, grad6, minimize = FALSE)
    if (lout$solution.type != "Optimal")
        stop("constraint set is empty (constraints are inconsistent)")
    if (qsign(lout$optimal.value) <= 0) {
        # step 2: discover which putative inequality constraints are
        # actually equality constraints
        if (stop.if.implied.equalities)
            stop("there are implied equalities, and you requested a stop")
        lin.time <- system.time(
            linout <- linearity(hrep1)
        )
        hrep1[linout, 1] <- "1"
    }

    stop.time.lin <- proc.time()

    hrep3 <- hrep1[hrep1[ , 1] == "1", , drop = FALSE]
    hrep4 <- hrep1[hrep1[ , 1] == "0", , drop = FALSE]

    # check to see that no variables are constrained equal to zero
    a3 <- hrep3[ , - c(1, 2), drop = FALSE]
    b3 <- hrep3[ , 2]
    solo <- apply(a3, 1, function(x) sum(x != "0") == 1)
    if (any(b3[solo] == "0")) {
        molo <- a3[solo, , drop = FALSE]
        nolo <- apply(molo, 1, function(x) seq(along = x)[x != "0"])
        colo <- paste(nolo, collapse = ", ")
        stop(paste("variable(s)", colo, "constrained to be exactly zero"))
    }

    # step 3: find V-representation of affine hull of constraint set
    # and use it to construct one-to-one affine transformation from
    # a vector space to the affine hull of constraint set
    # call this map from "new coordinates" (NC) to "old coordinates (OC)
    # also find constraint set in NC that maps one-to-one and onto
    # constraint set in OC using this one-to-one affine transformation

    basis.time <- system.time(
        vrep3 <- scdd(hrep3, representation = "H")$output
    )
    is.line <- vrep3[ , 1] == "1" & vrep3[ , 2] == "0"
    is.point <- vrep3[ , 1] == "0" & vrep3[ , 2] == "1"
    if (! all(is.point | is.line))
        stop("unexpected V-representation of affine hull of constraint set")
    if (sum(is.point) != 1)
        stop("unexpected V-representation of affine hull of constraint set")
    if (sum(is.line) == 0)
        stop("constraint set consists of single point or is empty")

    foo <- vrep3[ , - c(1, 2), drop = FALSE]
    origin <- foo[is.point, ]
    basis <- foo[is.line, , drop = FALSE]
    basis <- t(basis)

    # at this point
    #     fred <- function(x) origin + basis %*% x
    # maps from new coordinates (NC) onto the affine hull
    #     of the constraint (a convex polytope) in original coordinated (OC)

    amat <- qneg(hrep4[ , - c(1, 2), drop = FALSE])
    bvec <- hrep4[ , 2]
    bvec <- qmq(bvec, qmatmult(amat, cbind(origin)))
    amat <- qmatmult(amat, basis)

    # at this point
    #     sally <- function(x) all(amat %*% x <= bvec)
    # is the indicator function of a convex polytope in NC that is
    # mapped one-to-one onto the constraint set in OC by the function
    # fred defined in the previous comment

    # step 4: find a point (to be initial point of Markov chain) that
    # is relative interior point of constraint set in NC

    start.time.rip <- proc.time()

    hrep5 <- cbind("0", bvec, qneg(amat), "-1")
    grad5 <- c(rep("0", ncol(amat)), "1")
    lout <- lpcdd(hrep5, grad5, minimize = FALSE)
    if (lout$solution.type != "Optimal" || qsign(lout$optimal.value) <= 0)
        stop("constraint set is empty (constraints are inconsistent)")
    x <- lout$primal.solution
    # rip is relative interior point of constraint set in NC
    rip <- x[- length(x)]

    stop.time.rip <- proc.time()

    origin <- q2d(origin)
    basis <- q2d(basis)
    bvec <- q2d(bvec)
    amat <- q2d(amat)
    rip <- q2d(rip)

    # rcdd ver 1.1-9 and above burn random numbers and change .Random.seed
    # (forced by R CMD check for R-3.2.0 and up)
    # so put .Random.seed back where it was at the beginning
    assign(".Random.seed", saveseed, .GlobalEnv)

    out <- hitrunHelper(alpha, rip, nbatch, blen, nspac,
        origin, basis, amat, bvec, outmat, debug)

    if (! missing(a1)) {
        out$a1 <- a1
        out$b1 <- b1
    }
    if (! missing(a2)) {
        out$a2 <- a2
        out$b2 <- b2
    }
    foo <- list(linearity = stop.time.lin - start.time.lin, basis = basis.time,
        relative.interior.point = stop.time.rip - start.time.rip,
        mcmc = out$time)
    bar <- foo[[1]]
    for (i in 2:length(foo)) bar <- bar + foo[[i]]
    out$split.time <- foo
    out$time <- bar
    class(out) <- "hitrun"
    return(out)
}

hitrunHelper <- function(alpha, initial, nbatch, blen, nspac,
        origin, basis, amat, bvec, outmat, debug)
{
    if (! exists(".Random.seed")) runif(1)
    saveseed <- .Random.seed

    stopifnot(is.numeric(alpha))
    stopifnot(is.finite(alpha))
    stopifnot(is.vector(alpha))
    stopifnot(length(alpha) >= 2)
    stopifnot(alpha > 0)

    stopifnot(is.numeric(initial))
    stopifnot(is.finite(initial))

    stopifnot(is.numeric(nbatch))
    stopifnot(length(nbatch) == 1)
    stopifnot(nbatch == round(nbatch))
    stopifnot(nbatch >= 1)
    stopifnot(is.numeric(blen))
    stopifnot(length(blen) == 1)
    stopifnot(blen == round(blen))
    stopifnot(blen >= 1)
    stopifnot(is.numeric(nspac))
    stopifnot(length(nspac) == 1)
    stopifnot(nspac == round(nspac))
    stopifnot(nspac >= 1)

    stopifnot(is.numeric(basis))
    stopifnot(is.finite(basis))
    stopifnot(is.matrix(basis))
    stopifnot(is.numeric(origin))
    stopifnot(is.finite(origin))
    stopifnot(is.vector(origin))
    stopifnot(nrow(basis) == length(alpha))
    stopifnot(ncol(basis) == length(initial))
    stopifnot(length(origin) == length(alpha))

    stopifnot(is.numeric(amat))
    stopifnot(is.finite(amat))
    stopifnot(is.matrix(amat))
    stopifnot(is.numeric(bvec))
    stopifnot(is.finite(bvec))
    stopifnot(is.vector(bvec))
    stopifnot(nrow(amat) == length(bvec))
    stopifnot(ncol(amat) == length(initial))

    if (! is.null(outmat)) {
        stopifnot(is.numeric(outmat))
        stopifnot(is.finite(outmat))
        stopifnot(is.matrix(outmat))
        stopifnot(ncol(outmat) == length(alpha))
    }

    stopifnot(is.logical(debug))
    stopifnot(length(debug) == 1)

    storage.mode(basis) <- "double"
    storage.mode(amat) <- "double"
    if (is.matrix(outmat))
        storage.mode(outmat) <- "double"
    out.time <- system.time(
    out <- .Call(C_hitrun, as.double(alpha), as.double(initial),
        as.integer(nbatch), as.integer(blen), as.integer(nspac),
        as.double(origin), basis, amat, as.double(bvec), outmat,
        as.logical(debug))
    )

    out$initial.seed <- saveseed
    out$final.seed <- .Random.seed
    out$time <- out.time
    out$alpha <- alpha
    out$nbatch <- nbatch
    out$blen <- blen
    out$nspac <- nspac
    out$origin <- origin
    out$basis <- basis
    out$amat <- amat
    out$bvec <- bvec
    out$outmat <- outmat
    out$batch <- t(out$batch)
    out$debug <- debug
    if (debug) {
        out$current <- t(out$current)
        out$proposal <- t(out$proposal)
        out$z <- t(out$z)
    }
    class(out) <- "hitrun"
    return(out)
}

