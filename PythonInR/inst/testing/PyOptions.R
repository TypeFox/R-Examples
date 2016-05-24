#' # pyOptions
require(testthat)
require(PythonInR)
invisible(capture.output(pyConnect()))

## Options
expect_that(names(pyOptions()), 
	equals(c("numpyAlias", "useNumpy", "pandasAlias", "usePandas", "winPython364")))
expect_that(pyOptions("numpyAlias"), equals("numpy"))
expect_that(pyOptions("numpyAlias", "np"), equals("np"))
expect_that(pyOptions("numpyAlias"), equals("np"))
expect_that(pyOptions("numpyAlias", "numpy"), equals("numpy"))


