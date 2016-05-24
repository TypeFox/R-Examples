context("getNames")

test_that("verifies getNames", {
   e1 <- new.env(parent = emptyenv())
   e2 <- new.env(parent = e1)
   e1[["a"]] <- 1
   e1[["b"]] <- 2
   e2[["a"]] <- 3  # duplication
   e2[["c"]] <- 4
   actual <- DYM:::getNames("obj", envir=e2)
   expect_equal(actual, c("a", "b", "c"))
})

test_that("not contains special names", {
   e <- new.env(parent = emptyenv())
   e[["x"]] <- 1
   e[[".y"]] <- 2   # hidden
   e[["%z%"]] <- 3  # operator
   actual <- DYM:::getNames("obj", envir=e)
   expect_equal(actual, c("x"))
})

test_that("verifies getNames with mode=export", {
   e1 <- new.env(parent = emptyenv())
   e2 <- new.env(parent = e1)
   e1[["a"]] <- 1
   e1[["b"]] <- 2
   e2[["a"]] <- 3  # duplication
   e2[["c"]] <- 4
   actual <- DYM:::getNames("export", envir=e2)
   expect_equal(actual, c("a", "c"))
})
