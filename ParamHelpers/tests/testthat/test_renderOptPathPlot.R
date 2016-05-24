context("renderOptPathPlot")
test_that("renderOptPathPlot", {
  # Test 1D-1D + short names + alpha + title
  ps0 = makeParamSet(
    makeNumericParam("x")
  )
  op0 = makeOptPathDF(par.set = ps0, y.names = c("y"), minimize = TRUE)
  X = rnorm(14)
  dob = c(rep(0, 5), 1:2)
  for (i in 1:7) {
    addOptPathEl(op0, x = list(
      x = X[i * 2]),
      y = X[ i * 2 - 1],
      dob = dob[i])
  }
  pl = renderOptPathPlot(op0, iter = 0, alpha = FALSE)
  pl = sapply(c(0,2), FUN = function(x) {
    renderOptPathPlot(op0, iter = x, short.x.names = "e", short.y.names = "f")
  })
  pl = plotOptPath(op0, iters = 0:2, pause = FALSE, title = "Optimization Path")


  # Test 2D(mixed)-3D + marked + ggplot.theme + log + size
  ps1 = makeParamSet(
    makeNumericParam("x"),
    makeDiscreteParam("z", values = list("a", "b"))
  )
  op1 = makeOptPathDF(par.set = ps1, y.names = c("y1", "y2", "y3"), minimize = c(TRUE, FALSE, TRUE))
  X = rnorm(100)
  Z = rnorm(25)
  Z = ifelse(Z < 0, "a", "b")
  dob = c(rep(0, 5), 1:20)
  for (i in 1:25) {
    addOptPathEl(op1, x = list(x = abs(X[i * 4 - 3]), z = Z[i]),
      y = c(y1 = X[i * 4 - 2], y2 = abs(X[i * 4 - 1]), y3 = abs(X[i * 4])), dob = dob[i])
  }
  pl = renderOptPathPlot(op1, iter = 3)
  pl = sapply(c(0, 20), FUN = function(x) {
    renderOptPathPlot(op1, iter = x, marked = "best",
      ggplot.theme = ggplot2::theme(legend.position = "bottom"),
      size.points = 1, size.lines = 3)
  })
  pl = sapply(c(0, 20), FUN = function(x) {
    renderOptPathPlot(op1, iter = x, marked = c(4, 10, 18), x.over.time = list(c("x"), c("z")),
      y.over.time = list(c("y1"), c("y2", "y3")), log = c("x", "y2", "y3"))
  })


  # Test 1D(discrete)-2D + marked + limits + short names + rest variables
  ps2 = makeParamSet(
    makeDiscreteParam("x", values = list("a", "b"))
  )
  op2 = makeOptPathDF(par.set = ps2, y.names = c("y1", "y2"), minimize = c(TRUE, TRUE),
    include.extra = TRUE)
  X = rep(c("a", "b"), 4)
  Y = rnorm(14)
  dob = c(rep(0, 5), 1:2)
  for (i in 1:5) {
    addOptPathEl(op2, x = list(x = X[i]), y = c(y1 = Y[i], y2 = Y[7 + i]), dob = dob[i],
      extra = list(extra.var = NA))
  }
  for (i in 6:7) {
    addOptPathEl(op2, x = list(x = X[i]), y = c(y1 = Y[i], y2 = Y[7 + i]), dob = dob[i],
      extra = list(extra.var = i))
  }
  pl = renderOptPathPlot(op2, iter = 2, x.over.time = list(c("x"), c("extra.var")),
    short.rest.names = c("extra"))
  pl = sapply(c(0,2), FUN = function(x) {
    renderOptPathPlot(op2, iter = x, marked = c(3),
      xlim = list(YSpace = c(-10, 10)),
      ylim = list(YSpace = c(-10, 10), XSpace = c(0, 10)),
      short.x.names = "variable1", short.y.names = c("y", "z"))
  })


  # Test 3D(mixed)-1D + colours + missing values + limits + short names + scale +
  #                     impute.scale + impute.value
  ps3 = makeParamSet(
    makeNumericParam("x"),
    makeNumericParam("y", requires = quote(x < 0.7)),
    makeDiscreteParam("z", values = list("a", "b", "c"), requires = quote(x > -0.7))
  )
  op3 = makeOptPathDF(par.set = ps3, y.names = c("y1"), minimize = c(TRUE))
  X = rnorm(7)
  X2 = rnorm(7)
  Y = rnorm(7)
  X2[X >= 0.7] = NA
  X3 = rep(c("a", "b", "c"), 5)
  X3[X <= -0.7] = NA
  dob = c(rep(0, 5), 1:2)
  for (i in 1:7) {
    addOptPathEl(op3, x = list(x = X[i], y = X2[i], z = X3[i]),
      y = c(y1 = Y[i]), dob = dob[i])
  }

  pl = sapply(c(0,2), FUN = function(x) {
    renderOptPathPlot(op3, iter = x, scale = "globalminmax", impute.scale = 2,
      impute.value = "miss")
  })
  pl = sapply(c(0,2), FUN = function(x) {
    renderOptPathPlot(op3, iter = x, xlim = list(YSpace = c(-0.5, 0.5)),
      short.x.names = c("a", "b", "c"), colours = c("black", "yellow", "orange", "green"),
      scale = "globalminmax")
  })


  # Test subsetting
  pl = renderOptPathPlot(op1, iter = 2, subset.obs = c(1,2,6,7))
  pl = renderOptPathPlot(op1, iter = 0, subset.vars = 1)
  pl = renderOptPathPlot(op1, iter = 0, subset.targets = 2:3)
})
