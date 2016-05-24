context("find*Point(s)")

test_that("they refuse equality constraints", {
  constr <- simplexConstraints(4)
  expect_error(findInteriorPoint(constr))
  expect_error(findExtremePoints(constr))
  expect_error(findVertices(constr))
})

test_that("they detect unbounded polytopes", {
  constr <- list(
    constr = matrix(c(-1, 0, 0, 0, 0,
                      4, -1, 0, 0, 0,
                      0, 4, -1, 1, 0,
                      0, 0, 1, -1.333333,
                      -1, 0, 0, 0, 0, 4), nrow=5, ncol=5),
    rhs = rep(0, 5),
    dir = rep("<=", 5))

  expect_error(findInteriorPoint(constr), "No solution")
  expect_error(findExtremePoints(constr), "No solution")
  expect_error(findVertices(constr), "Failed to enumerate vertices")
})
