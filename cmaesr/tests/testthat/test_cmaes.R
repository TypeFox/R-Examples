context("CMA-ES run")

test_that("CMA-ES finds optimum of some BBOB functions", {
	max.iters = 100L
	lambda = 50L
	sigma = 1.5
  dims = c(2, 3, 4, 5, 10, 15, 20)

  # accepted tolerance value for parameter and fitness values
  tol = 0.05
  fun.generators = c(makeSphereFunction, makeAckleyFunction, makeDoubleSumFunction)
  stop.ons = c(list(stopOnMaxIters(max.iters)), getDefaultStoppingConditions())

  for (generator in fun.generators) {
    for (dim in dims) {
      fn = do.call(generator, list(dim))
      par.set = getParamSet(fn)
      opt = getGlobalOptimum(fn)
      lb = getLower(par.set)[1L]; ub = getUpper(par.set)[1L]

      res = cmaes(
        fn,
        start.point = runif(dim, min = lb, max = ub),
        monitor = NULL,
        control = list(lambda = lambda * dim, sigma = sigma, stop.ons = stop.ons)
      )

      expect_true(abs(res$best.fitness - opt$value) < tol,
        info = sprintf("Desired fitness level not reached for dim = %i and function '%s'", dim, getName(fn)))
      expect_true(sum((res$best.param - opt$param)^2) < tol,
        info = sprintf("Desired parameter approximation not reached for dim = %i and function '%s'", dim, getName(fn)))
    } # dims
  } # fun.generators
})

test_that("CMA-ES works on Sphere with default parameters", {
  # accepted tolerance value for parameter and fitness values
  tol = 0.1
  max.iters = 100L

  for (dim in c(2, 3, 5)) {
    fn = makeSphereFunction(dim)
    res = cmaes(
      fn,
      monitor = NULL,
      control = list(
        lambda = dim * 2 * 10,
        stop.ons = c(list(stopOnMaxIters(max.iters)), getDefaultStoppingConditions())
      )
    )
    expect_true(is.numeric(res$best.fitness))
    expect_true(all(is.numeric(res$best.param)))
    expect_true(res$best.fitness < tol, info = sprintf("For '%s' the desired fitness level was not reached.", getName(fn)))
  }
})

test_that("CMA-ES stops on invalid input", {
  control = list(sigma = 1)

  # multi-objective functions not supported
  fn = makeZDT1Function(2L)
  expect_error(cmaes(fn, monitor = NULL, control = control))

  # noisy functions not allowed
  fn = makeSphereFunction(2L)
  attr(fn, "noisy") = TRUE
  expect_error(cmaes(fn, monitor = NULL, control = control))

  # infinite bounds
  fn = makeSphereFunction(2L)
  attr(fn, "par.set") = makeNumericParamSet("x", len = 2L, lower = -Inf, upper = Inf)
  expect_error(cmaes(fn, monitor = NULL, control = control), regexp = "bounds")

  # negative weights
  fn = makeSphereFunction(2L)
  control2 = control
  control2$mu = 10
  control2$weights = runif(control2$mu)
  control2$weights[c(1, 3)] = -control2$weights[c(1, 3)]
  control2$stop.ons = getDefaultStoppingConditions()
  expect_error(cmaes(fn, monitor = NULL, control = control2), regexp = "negative")

  # missing stopping conditions
  fn = makeSphereFunction(2L)
  control = list(stop.ons = NULL, stop.ons = getDefaultStoppingConditions())
  expect_error(cmaes(fn, monitor = NULL, control = control), regexp = "stopping condition")

  # invalid "short name" as restart trigger
  control = list(restart.triggers = c("invalid_trigger"), stop.ons = getDefaultStoppingConditions())
  expect_error(cmaes(fn, monitor = NULL, control = control), regexp = "no stopping condition")

  # mixed functions
  fn = makeSingleObjectiveFunction(
    name = "Mixed",
    fn = function(x) {
      return(x$x^2 + as.numeric(x$y == "a"))
    },
    par.set = makeParamSet(
      makeNumericParam("x", lower = -10, upper = 10),
      makeDiscreteParam("y", values = c("a", "b"))
    )
  )
  expect_error(cmaes(fn, monitor = NULL, control = control))
})

test_that("CMA-ES computes reasonanable results on noiseless 2D BBOB test set", {
  # check all functions
  fids = 1:24
  dims = 2
  lambda = 200L
  tol = 0.5
  max.iters = 200L

  for (fid in fids) {
    for (dim in dims) {
      # skip the hardest (very multimodal) functions
      if (fid %in% c(24)) {
        next
      }
      lambda2 =  ifelse (fid %in% c(4, 5, 16, 23, 24), lambda * 4, lambda)
      fn = makeBBOBFunction(fid = fid, iid = 1L, dimension = dim)
      par.set = getParamSet(fn)
      opt = getGlobalOptimum(fn)
      lb = getLower(par.set)[1L]; ub = getUpper(par.set)[1L]
      control = list(
        sigma = (ub - lb) / 2,
        lambda = lambda2,
        stop.ons = c(list(stopOnMaxIters(max.iters)), getDefaultStoppingConditions())
      )

      res = cmaes(
        fn,
        control = control,
        monitor = NULL
      )
      expect_true(is.numeric(res$best.fitness))
      expect_true(abs(res$best.fitness - opt$value) < tol,
        info = sprintf("Desired fitness level not reached for dim = %i and function '%s'", dim, getName(fn)))
    }
  }
})

test_that("IPOP-CMA-ES works", {
  # pretty stupid example here to check if restarts are triggered:
  # We run CMA-ES for 5000 generations on Ackley and do trigger a restart
  # for all default stopping conditions.
  max.restarts = 3L
  max.iters = 5000L

  fn = makeAckleyFunction(2L)
  par.set = getParamSet(fn)
  lb = getLower(par.set); ub = getUpper(par.set)
  control = list(
    sigma = (ub[1L] - lb[1L]) / 2,
    lambda = 10L,
    stop.ons = c(list(stopOnMaxIters(max.iters)), getDefaultStoppingConditions()),
    max.restarts = 3L,
    restart.triggers = c("indefCovMat", "conditionCov", "noEffectAxis", "noEffectCoord", "tolX")
  )
  res = cmaes(fn, control = control, monitor = NULL)
  expect_true(is.numeric(res$best.fitness))
  expect_true(is.integer(res$n.restarts))
  expect_equal(res$n.restarts, max.restarts)
})

test_that("CMA-ES finds optimum even if it is located on the edge of the feasible region", {
  tol = 0.01

  # generate sphere function with slightly modified parameter set (optimum on the border)
  fn = makeSingleObjectiveFunction(
    name = "Evil Sphere",
    fn = function(x) {
      if (x[1L] < 0) -(sum(x)^2) else sum(x)^2
    },
    par.set =  makeNumericParamSet(
      len = 2L,
      id = "x",
      lower = c(0, -5.12),
      upper = c(5.12, 5.12),
      vector = TRUE
    ),
    global.opt.param = c(0, 0),
    global.opt.value = 0
  )

  res = cmaes(fn, monitor = NULL)
  expect_true(abs(res$best.fitness - getGlobalOptimum(fn)$value) < tol,
    info = sprintf("Did not find optimum on the ridge of feasible space"))
})
