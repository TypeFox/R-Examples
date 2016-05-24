#' # Basics
require(testthat)
require(PythonInR)
invisible(capture.output(pyConnect()))

#' ## pyDir
expect_that(intersect(pyDir(), "__name__"), equals("__name__"))
expect_that(intersect(pyDir("sys"), "version"), equals("version"))

#' ## pyHelp
expect_that(pyHelp("abs"), prints_text("built-in function"))

#' ## pyType
expect_that(pyType("dict()"), equals("dict"))
