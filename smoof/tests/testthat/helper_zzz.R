expectIsSmoofFunction = function(obj, generator) {
  expect_is(obj, c("smoof_function", "function"), info = "No snoof function generated '%s'.", generator)
}

expectGlobalOptimum = function(fun, generator) {
  op = getGlobalOptimum(fun)
  op.df = op$param
  for (i in 1:nrow(op.df)) {
    param = op.df[i, ]
    comp.op = fun(param)
    expect_true(abs(comp.op - op$value) < 0.01, info = sprintf("%i-th Global optimum does not correspond to given value
      for function '%s'! IS: %.4f, SHOULD BE: %.4f", i, generator, comp.op, op$value), label = generator)
  }
}

checkGGPlot = function(pl, title, xlab, ylab) {
  expect_is(pl, "gg")
  expect_is(pl, "ggplot")
  expect_equal(pl$labels$title, title)
  expect_equal(as.character(pl$labels$x), xlab)
  expect_equal(as.character(pl$labels$y), ylab)
}
