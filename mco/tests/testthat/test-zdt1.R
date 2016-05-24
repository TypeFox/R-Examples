## Check convergence towards the true PF which is given by
##
##  y_2 = 1 - sqrt(y_1)
context("zdt1")

set.seed(42)
r <- nsga2(zdt1, 30, 2,
           lower.bounds=rep(0, 30), upper.bounds=rep(1, 30),
           generations=1000)
y1 <- r$value[, 1]
y2 <- r$value[, 2]

test_that("Above Pareto front", {
  ## PF given by y_2 = 1 - sqrt(y_1)
  expect_true(all(y2 >= 1 - sqrt(y1)))
})

test_that("close to Pareto front", {
  ## Limit the maximum relative distance from the front
  rel_d_to_front <- max((y2 - (1 - sqrt(y1))) / max(y1, y2))
  expect_true(all(rel_d_to_front < 0.05))
})
