context("conversion of minimization to maximization problem and vice versa")

test_that("min<->max conversion works as expected", {
  # test function set
  funs = list(
    makeSphereFunction(1L),
    makeAckleyFunction(1L),
    makeRastriginFunction(1L)
  )

  checkSumOfMinMaxValuesVanishes = function(fun1, fun2, fun.name) {
    # generate some random parameters within the box constraints
    lb = getLowerBoxConstraints(fun)
    ub = getUpperBoxConstraints(fun)
    rnds = runif(50, min = lb, max = ub)

    # build sum (since f(x) + (-f(x)) = 0 should hold)
    sums = sapply(rnds, function(x) fun1(x) + fun2(x))
    expect_true(all(sums < 0.0001), info = sprintf("Function value differences are not equal to 0 for function %s.", getName(fun)))
    expect_true(shouldBeMinimized(fun1) != shouldBeMinimized(fun2))
  }

  for (fun in funs) {
    # convert to maximization problem and check if function values changed
    fun2 = convertToMaximization(fun)
    checkSumOfMinMaxValuesVanishes(fun, fun2, getName(fun))

    # convert back
    expect_error(convertToMaximization(fun2))
    fun3 = convertToMinimization(fun2)
    checkSumOfMinMaxValuesVanishes(fun2, fun3, getName(fun))
  }
})
