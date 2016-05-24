context("main functions")
test.sisal <- function() {
    X <- cbind(sine=sin((1:100)/5),
               linear=seq(from=-1, to=1, length.out=100),
               matrix(rnorm(200), 100, 2,
                      dimnames=list(NULL, paste("random", 1:2, sep="."))))
    y <- drop(X %*% c(3, 10, 1, 0) + rnorm(100))
    foo <- sisal(X, y, Mtimes = 50, kfold = 5, verbose = 0)
    test_that("Class of result is correct", {
        expect_true(inherits(foo, "sisal"))
    })
    sisal_names <-
        c("L.f", "L.v", "E.tr", "s.tr", "E.v", "L.f.nobranch", "L.v.nobranch",
          "E.tr.nobranch", "s.tr.nobranch", "E.v.nobranch", "n.evaluated",
          "edges", "vertices", "vertices.logical", "vertex.data", "var.names",
          "n", "d", "n.missing", "n.clean", "lm.L.f", "lm.L.v", "lm.full",
          "magic.L.f", "magic.L.v", "magic.full", "mean.y", "sd.y",
          "zeroRange.y", "mean.X", "sd.X", "zeroRange.X", "constant.X",
          "params", "pairwise.points", "pairwise.wins", "pairwise.preferences",
          "pairwise.rank", "path.length", "nested.path", "nested.rank",
          "branching.useful", "warnings", "n.warn")
    test_that("Result has the right items", {
        expect_named(foo, sisal_names)
    })
    ## TODO: more tests
}
test.sisal()
test.testSisal <- function() {
    ## TODO: real tests
    test_that("Dummy test passes", {
        expect_true(TRUE)
    })
}
test.testSisal()
test.bootMSE <- function() {
    ## TODO: real tests
    test_that("Dummy test passes", {
        expect_true(TRUE)
    })
}
test.bootMSE()
