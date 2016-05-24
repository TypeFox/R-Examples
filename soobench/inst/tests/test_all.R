check_fun <- function(fn) {
  context(function_id(fn))
  X <- random_parameters(100, fn)
  Y1 <- apply(X, 2, fn)
  Y2 <- fn(X)
  expect_equal(Y1, Y2)
  opt <- global_minimum(fn)
  if (is.list(opt$par)) {
    for (par in opt$par) {
      expect_that(fn(par),
                  equals(opt$value, sprintf("opt of %s",
                                            function_id(fn))))
    }
  } else {
    expect_that(fn(opt$par),
                equals(opt$value, sprintf("opt of %s", function_id(fn))))
  }
}

check_fun(goldstein_price_function())
check_fun(branin_function())

for (dim in c(2, 3, 5, 10, 20)) {
  check_fun(ackley_function(dim))
  check_fun(double_sum_function(dim))
  check_fun(ellipsoidal_function(dim))
  check_fun(griewank_function(dim))
  check_fun(kotancheck_function(dim))
  check_fun(mexican_hat_function(dim))
  check_fun(rastrigin_function(dim))
  check_fun(rosenbrock_function(dim))
  ##check_fun(schwefel_function(dim))
  check_fun(sphere_function(dim))
  check_fun(weierstrass_function(dim))
}

for (dim in c(2, 3, 5, 10, 20)) {
  for (fid in 1:24) {
    for (tid in 1:15) {
      check_fun(bbob2009_function(dim, fid, tid))
    }
  }
}
