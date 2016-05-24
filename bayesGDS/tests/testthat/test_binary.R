
context("binary")

test_that("small", {

    require(sparseHessianFD)
    require(sparseMVN)
    require(trustOptim)

    set.seed(123)

    rmvn.sparse.wrap <- function(n.draws, params) {
        rmvn.sparse(n.draws, params[["mean"]], params[["CH"]], prec=TRUE)
    }
    dmvn.sparse.wrap <- function(d, params) {
        dmvn.sparse(d, params[["mean"]], params[["CH"]], prec=TRUE)
    }



    data(binary_small)
    binary <- binary_small #rename for brevity

    N <- length(binary[["Y"]])
    k <- NROW(binary[["X"]])
    q <- k
    nvars <- as.integer(N*k + q)
    priors <- list(inv.Sigma = diag(k), ##rWishart(1,k+5,diag(k))[,,1],
                   inv.Omega = diag(k))

    start <- rnorm(nvars) ## random starting values
    f <- binary.f(start, data=binary, priors=priors)
    df <- binary.grad(start, data=binary, priors=priors)
    d2f <- binary.hess(start, data=binary, priors=priors)

    hs <- Matrix(0, nvars, nvars)
    for (i in 1:(N + 1)) {
        rng <- ((i-1)*k+1):(k*i)
        hs[rng, rng] <- tril(Matrix(1,k,k)) ## lower triangle
    }
    hs[N*k + 1:q, 1:(N*k)] <- 1 ## bottom margin
    hsNZ <- Matrix.to.Coord(hs)

    FD <- sparseHessianFD(start, binary.f, binary.grad,
                          rows=hsNZ[["rows"]], cols=hsNZ[["cols"]],
                          data=binary, priors=priors)

    ##----usingFD
    f <- FD$fn(start)
    df <- FD$gr(start)
    hess <- FD$hessian(start)

    expect_equivalent(hess, d2f)

    opt <- trust.optim(start, fn=FD$fn, gr = FD$gr, hs = FD$hessian,
                       method = "Sparse",
                       control = list(
                           start.trust.radius=5, stop.trust.radius = 1e-7,
                           prec=1e-7, report.precision=1,
                           maxit=500, preconditioner=1,
                           function.scale.factor=-1,
                           report.freq = 10000
                           )
                       )

    theta.star <- opt[["solution"]]
    hess <- opt[["hessian"]]

    ##----propParams
    scale <- .96
    chol.hess <- Cholesky(-scale*hess)
    prop.params <- list(mean = theta.star, CH = chol.hess)

    M <- 10000  ## proposal draws
    log.c1 <- FD$fn(theta.star)
    log.c2 <- dmvn.sparse.wrap(theta.star, prop.params)
    draws.m <- as(rmvn.sparse.wrap(M,prop.params),"matrix")
    log.post.m <- plyr::aaply(draws.m, 1, FD$fn)
    log.prop.m <- dmvn.sparse.wrap(draws.m, params=prop.params)
    log.phi <- log.post.m - log.prop.m + log.c2 - log.c1
    valid.scale <- all(log.phi <= 0)
    expect_true(valid.scale)


##----sampleGDS_serial
    n.draws <- 3  ## total number of draws needed
    max.tries <- 100000  ## to keep sample.GDS from running forever
    draws <- sample.GDS(n.draws = n.draws,
                        log.phi = log.phi,
                        post.mode = theta.star,
                        fn.dens.post = FD$fn,
                        fn.dens.prop = dmvn.sparse.wrap,
                        fn.draw.prop = rmvn.sparse.wrap,
                        prop.params = prop.params,
                        report.freq = 1, announce=TRUE)

    expect_true(all(draws$gt.1==0))
    expect_false(any(is.na(draws$counts)))

    LML <- get.LML(counts=draws$counts,
                   log.phi=log.phi,
                   post.mode=theta.star,
                   fn.dens.post= FD$fn,
                   fn.dens.prop=dmvn.sparse.wrap,
                   prop.params=prop.params)
    expect_true(is.finite(LML))
    expect_true(LML<0)

})
