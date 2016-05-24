
context("Cleaning up on error")

test_that("We clean up if a package cannot be installed", {

  lib_dir <- tempfile()
  try(
    silent = TRUE,
    suppressWarnings(
      pkgs <- make_packages(
        lib_dir = lib_dir,
        foo1 = { f <- function() print("hello!") ; d <- 1:10 },
        foo2 = { "syntax error" + 10 + "here" }
      )
    )
  )

  expect_false("package:foo1" %in% search())
  expect_false("package:foo2" %in% search())
  expect_false(file.exists(lib_dir))

})

test_that("We clean up if a package cannot be loaded", {

  lib_dir <- tempfile()
  try(
    silent = TRUE,
    suppressWarnings(
      pkgs <- make_packages(
        lib_dir = lib_dir,
        foo1 = { f <- function() print("hello!") ; d <- 1:10 },
        foo2 = { .onLoad <- function() { "syntax error" + 10 + "here" } }
      )
    )
  )

  expect_false("package:foo1" %in% search())
  expect_false("package:foo2" %in% search())
  expect_false(file.exists(lib_dir))

})
