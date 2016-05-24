context("test some getter and setter functions")

test_that("hasTags works as expected", {
  # check that hasTags works for smoof_functions
  fun.names = as.list(filterFunctionsByTags(c("multimodal")))
  funs = unlist(lapply(fun.names, function(fun.name) {
    makeFunctionsByName(fun.name, dimensions = 2L)
  }))

  for (fun in funs) {
    expect_true(hasTags(fun, "multimodal"), info = sprintf("Function '%s' has
      no tag multimodal, but it should have.", getName(fun)))
  }

  # check that hasTags works for smoof_generators
  expect_true(hasTags(makeAckleyFunction, "multimodal"))
  expect_false(hasTags(makeSphereFunction, "multimodal"))

  # check that hasTags works for strings
  expect_true(hasTags("Ackley", "multimodal"))
  expect_false(hasTags("Sphere", "multimodal"))
})


test_that("getLocalOptimum and getGlobalOptimum work as expected on MPM2", {
  # check that hasTags works for smoof_functions
  n.peaks = c(1, 5, 30)
  dims = c(2, 3, 10)
  for (d in dims) {
    funs = lapply(n.peaks, function(peaks) makeMPM2Function(n.peaks = peaks,
      dimensions = d, topology = "funnel", seed = 1, rotated = TRUE, peak.shape = "ellipse"))

    for (i in seq_along(funs)) {
      fun = funs[[i]]
      local.opt = getLocalOptimum(fun)
      global.opt = getGlobalOptimum(fun)
      # check that each peak has a parameter and objective value
      expect_equal(length(local.opt$params), n.peaks[i])
      expect_equal(length(local.opt$values), n.peaks[i])
      # length of x-values of local optima should equal the problem dimension
      expect_equal(unique(vapply(local.opt$params, length, integer(1L))), d)
      expect_equal(unique(vapply(local.opt$values, length, integer(1L))), 1)
      # check that global optimum is the best of the local optima
      best = ifelse(global.opt$is.minimum, min, max)
      indices = which(best(unlist(local.opt$values)) == global.opt$value)
      for (ind in indices) {
        expect_equal(unlist(local.opt$values)[ind], global.opt$value)
      }
      expect_true(any(vapply(indices, function(ind)
        all.equal(local.opt$params[[ind]], global.opt$param), logical(1L))))
    }
  }
})
