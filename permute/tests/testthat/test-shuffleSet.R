library("testthat")
library("permute")

context("Testing shuffleSet()")

## test that shuffleSet interleves the separate block-level
## permutations correctly back into the original ordering
## This always generates odd, even, odd, ..., numbered vector
## of observations, hence when we take the modulus 2 we get
## a vector of 1,0,1,0,1,...
test_that("shuffleSet interleves block-level perms correctly", {
    gr <- factor(rep(1:2, length=20))
    ctrl <- how(nperm = 5, blocks = gr)
    p <- shuffleSet(20, control = ctrl) %% 2
    y <- rep(c(1L, 0L), length.out = ncol(p))
    nc <- ncol(p)
    for (i in seq_len(nrow(p))) {
        expect_that(p[i, ], equals(y))
    }
})

## test that nset permutations are always returned if
## make = FALSE in how()
test_that( "shuffleSet returns exactly nset permutations when make == FALSE", {
    ## simple random permutation
    h <- how(make = FALSE)
    ss <- shuffleSet(n = 4, nset = 10, control = h)
    expect_that(nrow(ss), equals(10))

    ## time series
    h <- how(within = Within(type = "series"), make = FALSE)
    ss <- shuffleSet(n = 20, nset = 15, control = h)
    expect_that(nrow(ss), equals(15))

})

## test that shuffleSet always returns a matrix, even for nset == 1
test_that("shuffleSet returns a matrix even for nset == 1", {
    h <- how()
    ss <- shuffleSet(25, nset = 1, control = h)
    expect_that(ss, is_a("matrix"))
})

test_that("shuffle can permute both plots and within in presence of blocks", {
    ## Example from @LindsayVass on github #9
    control <- how(within = Within(type = "free"),
                   plots = Plots(strata = rep(gl(2,7),2), type = "free"),
                   blocks = gl(2, 14))
    permSet <- shuffleSet(28, 100, control = control)
    expect_that(nrow(permSet), is_identical_to(100L))
    expect_that(ncol(permSet), is_identical_to(28L))
    expect_that(permSet, is_a("permutationMatrix"))
    expect_that(permSet, is_a("matrix"))
})

test_that("print method for permutationMatrix works", {
    h <- how()
    perms <- shuffleSet(10, nset = 10, control = h)
    expect_output(print(perms), regexp = "No. of Permutations:")

    h <- how(blocks = gl(5,10))
    perms <- shuffleSet(50, nset = 20, control = h)
    expect_output(print(perms), regexp = "Restricted by Blocks:")
    expect_output(print(perms), regexp = "blocks;")

    h <- how(plots = Plots(strata = gl(5,10)))
    perms <- shuffleSet(50, nset = 20, control = h)
    expect_output(print(perms), regexp = "Restricted by Plots:")
    expect_output(print(perms), regexp = "plots;")

    h <- how(plots = Plots(strata = gl(10,10)),
             blocks = gl(2, 50))
    perms <- shuffleSet(100, nset = 20, control = h)
    expect_output(print(perms), regexp = "Restricted by Plots:")
    expect_output(print(perms), regexp = "plots & blocks;")

    h <- how(within = Within(type = "series", mirror = TRUE))
    perms <- shuffleSet(10, nset = 20, control = h)
    expect_output(print(perms), regexp = "; mirrored")

    h <- how(within = Within(type = "series", constant = TRUE),
             plots = Plots(strata = gl(2, 5), type = "series", mirror = TRUE))
    perms <- shuffleSet(10, nset = 20, control = h)
    expect_output(print(perms), regexp = "; same permutation")
})
