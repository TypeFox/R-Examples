context("plapply()")

test_that("plapply() results are reproducible", {
    if (.Platform$OS.type == "windows")
        skip("causes R CMD check to hang on Windows (but works interactively)")
    res1 <- plapply(c(1, 1), rnorm, .parallel = 2, .seed = 1, .verbose = FALSE)
    res2 <- plapply(c(1, 1), rnorm, .parallel = 2, .seed = 1, .verbose = FALSE)
    expect_that(res1, is_identical_to(res2))
})
