#' # pyAttach
require(testthat)
require(PythonInR)
invisible(capture.output(pyConnect()))

expect_that(pyExec("import os"), equals(0))

## attach to global
## ----------------
## attach the function getcwd from the module os to R.
expect_that(pyAttach("os.getcwd", .GlobalEnv), is_a("environment"))
expect_that(os.getcwd(), is_a("character"))
## attach the string object os.name to R
expect_that(pyAttach("os.name", .GlobalEnv), is_a("environment"))
expect_that(pyGet("os.name"), equals(os.name))

## Since os.name is attached to the globalenv it can be set without using
## the global assignment operator
os.name <<- "Hello Python from R!"
expect_that(pyGet("os.name"), equals(os.name))
