context("assignment")

test_that("assignment", {
  expect_null({
    declare(x,y,z)
    numeric(x) <- rnorm(5)
    integer(y) <- 2L
    logical(z) <- TRUE
    check(x, class="numeric") <- rnorm(5)
    NULL
  })
  expect_error({
    declare(x=numeric())
    integer(x) <- 10L
  })
  expect_error({
    declare(x=numeric())
    numeric(x) <- TRUE
  })
  expect_error({
    declare(x=numeric())
    numeric(x, length = 10) <- 1:9
  })
  expect_error({
    declare(x=numeric())
    cond <- function(x) mean(x) <= 5
    numeric(x, cond) <- 1:10
  })
  expect_error({
    declare(x=numeric())
    cond <- function(m) {
      function(x) mean(x) <= m
    }
    numeric(x, cond(5)) <- 1:10
  })
})
