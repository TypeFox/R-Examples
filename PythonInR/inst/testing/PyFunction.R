#' # pyImport
require(testthat)
require(PythonInR)
invisible(capture.output(pyConnect()))

expect_that(pyExec("import os"), equals(0))
fun <- pyFunction("os.getcwd")
expect_that(class(fun), equals("pyFunction"))
expect_that(fun(), equals(getwd()))
