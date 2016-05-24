context("Seed point generation tests")

test_that("seed point is inside the polytope", {
  # 1 > 4 > 5 > 6 > {2, 3}
  constr <- structure(list(constr = structure(c(-1, 0, 0, 0, 0, 0, 0, 0,
    1, 0, 0, 0, 0, 0, 1, 1, -1, 0, 0, 0, 0, 1, -1, 0, 0, 0, 0, 1,
    -1, -1), .Dim = 5:6), rhs = c(0, 0, 0, 0, 0), dir = c("<=", "<=",
    "<=", "<=", "<=")), .Names = c("constr", "rhs", "dir"))
  n <- ncol(constr$constr)

  transform <- simplex.createTransform(n)
  constr <- simplex.createConstraints(transform, constr)
  seedPoint <- createSeedPoint(constr, homogeneous=TRUE)

  expect_true(all((constr$constr %*% seedPoint) <= constr$rhs))
})

test_that("randomized seed point generation works", {
  n <- 4
  ord <- c(4,2,1,3)
  pairs <- cbind(ord[1:(n-1)], ord[2:n])
  constr <- mergeConstraints(
    apply(pairs, 1, function(pair) {
      ordinalConstraint(n, pair[1], pair[2])
    }))
  transform <- simplex.createTransform(n)
  constr <- simplex.createConstraints(transform, constr)

  seedPoint <- createSeedPoint(constr, homogeneous=TRUE)
  expect_true(all((constr$constr %*% seedPoint) <= constr$rhs))

  s1 <- createSeedPoint(constr, method="slacklp", randomize=TRUE, homogeneous=TRUE)
  expect_true(all((constr$constr %*% s1) <= constr$rhs))

  s2 <- createSeedPoint(constr, method="slacklp", randomize=TRUE, homogeneous=TRUE)
  expect_true(all((constr$constr %*% s2) <= constr$rhs))
  expect_false(isTRUE(all.equal(s1, s2)))
  expect_false(isTRUE(all.equal(seedPoint, s2)))
})
