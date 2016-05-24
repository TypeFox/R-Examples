context("plotEAF")

test_that("plotEAF works", {

  # Build 3 artifical opt.paths with 2 crits
  makeArtOP = function() {
    ps = makeNumericParamSet(len = 1L)
    op = makeOptPathDF(par.set = ps, y.names = c("y1", "y2"),
      minimize = c(TRUE, TRUE))
    X = rnorm(60)
    dob = c(rep(0, 5), 1:15)
    for (i in 1:20) {
      addOptPathEl(op, x = list(
        x = X[i * 3]),
        y = X[ i * 3 - 1:2],
        dob = dob[i])
    }
    op
  }
  # build aritificial list of opt pathes
  opt.paths = list(
    algo1 = list(makeArtOP(), makeArtOP(), makeArtOP()),
    algo2 = list(makeArtOP(), makeArtOP()),
    algo3  = list(makeArtOP()))

  # plot eaf and check returned data frame
  res = plotEAF(opt.paths)

  algo.names = c("algo1", "algo2", "algo3")

  expect_true(is.data.frame(res))
  expect_true(setequal(colnames(res), c("y1", "y2", ".algo", ".repl")))
  expect_true(all(res$.algo %in% algo.names))
  expect_true(all(res$.repl >= 1) && all(res$.repl <= 3))
})
