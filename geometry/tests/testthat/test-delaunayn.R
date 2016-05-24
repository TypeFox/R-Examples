context("delaunayn")
test_that("delaunayn produces the correct output", {
  ## Create points that, when passed to Qhull with the Qt option,
  ## would give degenerate simplices - thanks to Bill Denney for
  ## example
  ps <- as.matrix(rbind(data.frame(a=0, b=0, d=0),
                        merge(merge(data.frame(a=c(-1, 1)),
                                    data.frame(b=c(-1, 1))),
                              data.frame(d=c(-1, 1)))))
  ts <- delaunayn(ps)
  expect_that(ts, is_a("matrix"))
  
  ## With full output, there should be a trinagulation, areas and
  ## neighbours and the sum of the ares should be 8
  ts.full <- delaunayn(ps, full=TRUE)
  expect_that(ts, equals(ts.full$tri))
  expect_that(length(ts.full$areas), equals(nrow(ts.full$tri)))
  expect_that(length(ts.full$neighbours), equals(nrow(ts.full$tri)))
  expect_that(sum(ts.full$area), equals(8))
  
  ## tsearchn shouldn't return a "degnerate simplex" error. 
  expect_that(tsearchn(ps, ts, cbind(1, 2, 4)), not(gives_warning("Degenerate simplices")))

  ## If the input matrix contains NAs, delaunayn should return an error
  ps <- rbind(ps, NA)
  expect_error(delaunayn(ps))

})

test_that("In the case of just one triangle, delaunayn returns a matrix", {
  pc  <- rbind(c(0, 0), c(0, 1), c(1, 0))
  pct <- delaunayn(pc)
  expect_that(pct, is_a("matrix"))
  expect_that(nrow(pct), equals(1))
  pct.full <- delaunayn(pc, full=TRUE)
  expect_that(pct.full$areas, equals(0.5))
})

test_that("In the case of a degenerate triangle, delaunayn returns a matrix with zero rows", {
  pc  <- rbind(c(0, 0), c(0, 1), c(0, 2))
  pct <- delaunayn(pc)
  expect_that(pct, is_a("matrix"))
  expect_that(nrow(pct), equals(0))
  pct.full <- delaunayn(pc, full=TRUE)
  expect_that(length(pct.full$areas), equals(0))
  expect_that(length(pct.full$neighbours), equals(0))
})

test_that("In the case of just one tetrahaedron, delaunayn returns a matrix", {
  pc  <- rbind(c(0, 0, 0), c(0, 1, 0), c(1, 0, 0), c(0, 0, 1))
  pct <- delaunayn(pc)
  expect_that(pct, is_a("matrix"))
  expect_that(nrow(pct), equals(1))
  pct.full <- delaunayn(pc, full=TRUE)
   expect_that(pct.full$areas, equals(1/6))
})

test_that("Output to file works", {
  ps <-  matrix(rnorm(3000), ncol=3)
  ps <-  sqrt(3)*ps/drop(sqrt((ps^2) %*% rep(1, 3)))
  pst <- delaunayn(ps, "QJ TO 'test1.txt'")
  expect_true(file.exists("test1.txt"))
})

test_that("The QJ option can give degenerate simplices", {
  ## Create degenerate simplex - thanks to Bill Denney for example
  ps <- as.matrix(rbind(data.frame(a=0, b=0, d=0),
                        merge(merge(data.frame(a=c(-1, 1)),
                                    data.frame(b=c(-1, 1))),
                              data.frame(d=c(-1, 1)))))

  ## The QJ option leads to on simplex being very small
  ts <- delaunayn(ps, "QJ")
  expect_warning(tsearchn(ps, ts, cbind(1, 2, 4)))
})


