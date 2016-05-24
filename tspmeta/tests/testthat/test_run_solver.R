context("run_solver")

if (interactive()) {
test_that("run_solver", {
  x = random_instance(size = 100)
  for (method in get_solvers()) {
    res = run_solver(x, method = method)
  }
})
}