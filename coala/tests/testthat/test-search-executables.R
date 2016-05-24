context("Searching Executables")

test_that("providing exes via ENVIR vars works", {
  exe <- tempfile("test_exe")
  Sys.setenv("COALATESTEXE" = exe)
  cat("test\n", file = exe)
  expect_equal(search_executable("test_exe", envir_var = "COALATESTEXE"), exe)
  unlink(exe)
  Sys.unsetenv("COALATESTEXE")
})


test_that("a warning is thrown if a non-exisiting exe is given via ENVIR var", {
  exe <- tempfile("test_exe")
  Sys.setenv("COALATESTEXE" = exe)
  expect_warning(search_executable("test_exe", envir_var = "COALATESTEXE"))
  Sys.unsetenv("COALATESTEXE")
})


test_that("finding executables in the current dir works", {
  current_wd <- getwd()
  setwd(tempdir())
  exe <- "coala_test_exe"
  cat("test\n", file = exe)
  expect_equal(search_executable("coala_test_exe"), file.path(getwd(), exe))
  unlink(exe)
  setwd(current_wd)
})
