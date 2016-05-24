library("testthat")
library("permute")

context("Testing allPerms()")

test_that("allPerms - blocks - within block free", {
    ## example data from Joris Meys from
    ## http://stackoverflow.com/a/21313632/429846
    thedata <- data.frame(score = c(replicate(4, sample(1:3))),
                          judge = rep(1:4, each = 3),
                          wine = rep.int(1:3, 4))

    ## without the observed permutation included
    hh <- how(within = Within("free"),
              blocks = factor(thedata$judge),
              complete = TRUE, maxperm = 1e9)
    nr <- nrow(thedata)
    np <- numPerms(nr, hh)
    p <- allPerms(nr, control = hh)
    expect_that(nrow(p), equals(np - 1)) ## default is to drop observed

    ## check no duplicate indices within rows
    dup <- any(apply(p, 1, function(x) any(duplicated(x))))
    expect_false(dup, info = "Blocks: even; within: free; no observed")

    ## with the observed permutation included
    hh <- how(within = Within("free"),
              blocks = factor(thedata$judge),
              complete = TRUE, maxperm = 1e9,
              observed = TRUE)
    p <- allPerms(nr, control = hh)
    expect_that(nrow(p), equals(np)) ## now includes observed

    ## check no duplicate indices within rows
    dup <- any(apply(p, 1, function(x) any(duplicated(x))))
    expect_false(dup, info = "Blocks: even; within: free; observed")
})

test_that("allPerms; blocks: within; block free - uneven block sizes", {
    fac <- factor(rep(1:3, times = c(2,2,4)))

    ## without the observed permutation included
    hh <- how(within = Within("free"),
              blocks = fac,
              complete = TRUE, maxperm = 1e9)
    ll <- length(fac)
    np <- numPerms(ll, hh)
    expect_that(np, equals(prod(factorial(2), factorial(2), factorial(4))))
    p <- allPerms(ll, control = hh)
    expect_that(nrow(p), equals(np - 1)) ## default is to drop observed

    ## check no duplicate indices within rows
    dup <- any(apply(p, 1, function(x) any(duplicated(x))))
    expect_false(dup, info = "Blocks: uneven; within: free; no observed")

    ## with the observed permutation included
    hh <- how(within = Within("free"),
              blocks = fac,
              complete = TRUE, maxperm = 1e9,
              observed = TRUE)
    p <- allPerms(ll, control = hh)
    expect_that(nrow(p), equals(np)) ## now includes observed

    ## check no duplicate indices within rows
    dup <- any(apply(p, 1, function(x) any(duplicated(x))))
    expect_false(dup, info = "Blocks: uneven; within: free; observed")
})

## testing plot-level permutations ------------------------------------
test_that("allPerms: plots; within: free; even: yes;", {
    fac <- rep(1:3, each = 3)

    hh <- how(plots = Plots(strata = fac),
              complete = TRUE, maxperm = 1e9)
    ll <- length(fac)
    np <- numPerms(ll, hh)
    p <- allPerms(ll, control = hh)
    expect_that(nrow(p), equals(np - 1), ## default is to drop observed
                info = "Check n all perms == numPerms output.")

    ## check no duplicate indices within rows
    dup <- any(apply(p, 1, function(x) any(duplicated(x))))
    expect_false(dup,
                 info = "Unique? Plots: even; within: free; no observed")

    ## with the observed permutation included
    hh <- how(within = Within("free"),
              plot = Plots(strata = fac),
              complete = TRUE, maxperm = 1e9,
              observed = TRUE)
    p <- allPerms(ll, control = hh)
    expect_that(nrow(p), equals(np)) ## now includes observed

    ## check no duplicate indices within rows
    dup <- any(apply(p, 1, function(x) any(duplicated(x))))
    expect_false(dup, info = "Unique? Plots: even; within: free; inc observed")
})

test_that("allPerms; plots: within; plot free - uneven plot sizes", {
    fac <- factor(rep(1:3, times = c(2,2,4)))

    ## without the observed permutation included
    hh <- how(within = Within("free"),
              plots = Plots(strata = fac),
              complete = TRUE, maxperm = 1e9)
    ll <- length(fac)
    np <- numPerms(ll, hh)
    expect_that(np, equals(prod(factorial(2), factorial(2), factorial(4))))
    p <- allPerms(ll, control = hh)
    expect_that(nrow(p), equals(np - 1)) ## default is to drop observed

    ## check no duplicate indices within rows
    dup <- any(apply(p, 1, function(x) any(duplicated(x))))
    expect_false(dup, info = "Plots: uneven; within: free; no observed")

    ## with the observed permutation included
    hh <- how(within = Within("free"),
              plots = Plots(strata = fac),
              complete = TRUE, maxperm = 1e9,
              observed = TRUE)
    p <- allPerms(ll, control = hh)
    expect_that(nrow(p), equals(np)) ## now includes observed

    ## check no duplicate indices within rows
    dup <- any(apply(p, 1, function(x) any(duplicated(x))))
    expect_false(dup, info = "Plots: uneven; within: free; observed")
})

test_that("allPerms; permuting plots only -- non-contiguous plots", {
    transect <- rep(gl(2,2), 2)
    ll <- length(transect)
    ctrl <- how(Within(type = "none"), Plots(type = "free", strata = transect))

    ## without observed
    ref <- matrix(c(3L,4L,1L,2L,7L,8L,5L,6L), nrow = 1, byrow = TRUE)
    perm <- allPerms(ll, ctrl)
    attr(perm, "control") <- NULL
    attr(perm, "observed") <- NULL
    class(perm) <- "matrix"
    expect_that(numPerms(ll, control = ctrl), equals(2L),
                info = "Number of permutations is wrong")
    expect_that(nrow(perm), equals(1L),
                info = "Number of rows in permutation matrix != 1")
    expect_identical(perm, ref)

    ## with observed
    setObserved(ctrl) <- TRUE
    ref <- matrix(c(1L,2L,3L,4L,5L,6L,7L,8L,
                    3L,4L,1L,2L,7L,8L,5L,6L), nrow = 2, byrow = TRUE)
    perm <- allPerms(ll, ctrl)
    perm <- as.matrix(perm)
    expect_that(numPerms(ll, control = ctrl), equals(2L),
                info = "Number of permutations is wrong")
    expect_that(nrow(perm), equals(2L),
                info = "Number of rows in permutation matrix != 2")
    expect_identical(perm, ref,
                     info = "All permutations doesn't match reference")
})

## Grid permutations
test_that("Can generate permutations from a grid design", {
    ## spatial grids within each level of plot, 3 x (4r x 4c)
    nr <- 4
    nc <- 4
    np <- 3 ## number of plots
    plots <- Plots(gl(np, prod(nr, nc)))
    CTRL <- how(plots = plots,
                within = Within(type = "grid", ncol = nc, nrow = nr))
    perms <- allPerms(prod(nr, nc, np), control = CTRL)
    nperms <- numPerms(prod(nr, nc, np), control = CTRL)

    expect_is(perms, "allPerms")
    expect_is(perms, "matrix")
    expect_equal(nperms, nrow(perms) + 1L)

    ## mirroring
    nr <- 3
    nc <- 3
    np <- 2 ## number of plots
    plots <- Plots(gl(np, prod(nr, nc)))
    CTRL <- how(plots = plots,
                within = Within(type = "grid", ncol = nc, nrow = nr,
                                mirror = TRUE))
    perms <- allPerms(prod(nr, nc, np), control = CTRL)
    nperms <- numPerms(prod(nr, nc, np), control = CTRL)

    expect_is(perms, "allPerms")
    expect_is(perms, "matrix")
    expect_equal(nperms, nrow(perms) + 1L)
})

test_that("grids with 2 columns only work correctly", {
    ## spatial grids within each level of plot, (4r x 2c)
    nr <- 4
    nc <- 2
    CTRL <- how(within = Within(type = "grid", ncol = nc, nrow = nr))
    perms <- allPerms(prod(nr, nc), control = CTRL)
    nperms <- numPerms(prod(nr, nc), control = CTRL)

    expect_is(perms, "allPerms")
    expect_is(perms, "matrix")
    expect_equal(nperms, nrow(perms) + 1L)

    ## spatial grids within each level of plot, 3 x (4r x 2c)
    nr <- 4
    nc <- 2
    np <- 3 ## number of plots
    plots <- Plots(gl(np, prod(nr, nc)))
    CTRL <- how(plots = plots,
                within = Within(type = "grid", ncol = nc, nrow = nr))
    perms <- allPerms(prod(nr, nc, np), control = CTRL)
    nperms <- numPerms(prod(nr, nc, np), control = CTRL)

    expect_is(perms, "allPerms")
    expect_is(perms, "matrix")
    expect_equal(nperms, nrow(perms) + 1L)
})


test_that("grids with mirroring & 2 columns only work correctly", {
    ## spatial grids within each level of plot, (4r x 2c)
    nr <- 4
    nc <- 2
    CTRL <- how(within = Within(type = "grid", ncol = nc, nrow = nr,
                                mirror = TRUE))
    perms <- allPerms(prod(nr, nc), control = CTRL)
    nperms <- numPerms(prod(nr, nc), control = CTRL)

    expect_is(perms, "allPerms")
    expect_is(perms, "matrix")
    expect_equal(nperms, nrow(perms) + 1L)

    ## spatial grids within each level of plot, 3 x (4r x 2c)
    nr <- 4
    nc <- 2
    np <- 3 ## number of plots
    plots <- Plots(gl(np, prod(nr, nc)))
    CTRL <- how(plots = plots,
                within = Within(type = "grid", ncol = nc, nrow = nr,
                                mirror = TRUE))
    perms <- allPerms(prod(nr, nc, np), control = CTRL)
    nperms <- numPerms(prod(nr, nc, np), control = CTRL)

    expect_is(perms, "allPerms")
    expect_is(perms, "matrix")
    expect_equal(nperms, nrow(perms) + 1L)
})

test_that("same grid permutation within plots", {
    ## spatial grids within each level of plot, 3 x (4r x 2c)
    nr <- 4
    nc <- 2
    np <- 3 ## number of plots
    plots <- Plots(gl(np, prod(nr, nc)))
    CTRL <- how(plots = plots,
                within = Within(type = "grid", ncol = nc, nrow = nr,
                                constant = TRUE))
    perms <- allPerms(prod(nr, nc, np), control = CTRL)
    nperms <- numPerms(prod(nr, nc, np), control = CTRL)

    expect_is(perms, "allPerms")
    expect_is(perms, "matrix")
    expect_equal(nperms, nrow(perms) + 1L)
})

test_that("same grid permutation within plots & mirroring", {
    ## spatial grids within each level of plot, 3 x (4r x 2c)
    nr <- 4
    nc <- 2
    np <- 3 ## number of plots
    plots <- Plots(gl(np, prod(nr, nc)))
    CTRL <- how(plots = plots,
                within = Within(type = "grid", ncol = nc, nrow = nr,
                                constant = TRUE, mirror = TRUE))
    perms <- allPerms(prod(nr, nc, np), control = CTRL)
    nperms <- numPerms(prod(nr, nc, np), control = CTRL)

    expect_is(perms, "allPerms")
    expect_is(perms, "matrix")
    expect_equal(nperms, nrow(perms) + 1L)
})

test_that("allPerms works with complex, but small, design", {
    h <- how(within = Within(type = "series", constant = TRUE),
             plots = Plots(strata = gl(2, 5), type = "series", mirror = TRUE))
    ap <- allPerms(10, control = h)
    expect_is(ap, "matrix")
    expect_equal(nrow(ap), 10 - 1L)
})
