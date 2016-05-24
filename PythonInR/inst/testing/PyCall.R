#' # pyCall
require(testthat)
require(PythonInR)
invisible(capture.output(pyConnect()))

expect_that(pyExec("import os"), equals(0))

#' ## get/set current working directory
## by default
expect_that(pyCall("os.chdir", args=list(getwd())), equals(NULL))
expect_that(gsub("\\", "/", pyCall("os.getcwd"), fixed=TRUE), equals(getwd()))

#' ## test builtins functions (__builtins__)
pyImport("sys")
expect_that(pyCall("abs", args=list(-5)), equals(5))
expect_that(pyCall("sum", args=list(1:5)), equals(sum(1:5)))
#' NOTE: since all integer variables are translated to long by default
#' the following code produces an error
expect_that(
    expect_that(
        pyCall("hex", args=list(255), namespace=builtinNsp),
        prints_text("TypeError")
    ),
    throws_error()
)
