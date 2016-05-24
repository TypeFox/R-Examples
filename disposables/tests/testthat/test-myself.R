
context("Disposable packages")

test_that("We can create and load them", {

  on.exit(dispose_packages(pkgs))
  pkgs <- make_packages(
    foo1 = { f <- function() print("hello!") ; d <- 1:10 },
    foo2 = { f <- function() print("hello again!") ; d <- 11:20 }
  )

  expect_output(foo1::f(), "hello!", fixed = TRUE)
  expect_output(foo2::f(), "hello again!", fixed = TRUE)
  expect_equal(foo1::d, 1:10)
  expect_equal(foo2::d, 11:20)

})

test_that("We can dispose packages", {

  on.exit(dispose_packages(pkgs))
  pkgs <- make_packages(
    foo1 = { f <- function() print("hello!") ; d <- 1:10 },
    foo2 = { f <- function() print("hello again!") ; d <- 11:20 }
  )

  ## Do nothing
  dispose_packages(pkgs, unattach = FALSE, delete = FALSE)
  expect_true("package:foo1" %in% search())
  expect_true("package:foo2" %in% search())
  expect_true("foo1" %in% loadedNamespaces())
  expect_true("foo2" %in% loadedNamespaces())
  expect_true("foo1" %in% dir(pkgs$lib_dir))
  expect_true("foo2" %in% dir(pkgs$lib_dir))

  ## Unattach
  dispose_packages(pkgs, unload = FALSE, delete = FALSE)
  expect_false("package:foo1" %in% search())
  expect_false("package:foo2" %in% search())
  expect_true("foo1" %in% loadedNamespaces())
  expect_true("foo2" %in% loadedNamespaces())
  expect_true("foo1" %in% dir(pkgs$lib_dir))
  expect_true("foo2" %in% dir(pkgs$lib_dir))

  ## Unload
  dispose_packages(pkgs, delete = FALSE)
  expect_false("package:foo1" %in% search())
  expect_false("package:foo2" %in% search())
  expect_false("foo1" %in% loadedNamespaces())
  expect_false("foo2" %in% loadedNamespaces())
  expect_true("foo1" %in% dir(pkgs$lib_dir))
  expect_true("foo2" %in% dir(pkgs$lib_dir))

  ## Delete
  dispose_packages(pkgs, delete_lib_dir = FALSE)
  expect_false("foo1" %in% dir(pkgs$lib_dir))
  expect_false("foo2" %in% dir(pkgs$lib_dir))
  expect_true(file.exists(pkgs$lib_dir))

  ## Clean up completely
  dispose_packages(pkgs)
  expect_false(file.exists(pkgs$lib_dir))

})

test_that("We unload a package if it is already loaded", {

  on.exit(dispose_packages(pkgs), add = TRUE)
  pkgs <- make_packages(
    foo1 = { f <- function() print("hello!") ; d <- 1:10 },
    foo2 = { f <- function() print("hello again!") ; d <- 11:20 }
  )

  on.exit(dispose_packages(pkgs2), add = TRUE)
  pkgs2 <- make_packages(
    foo1 = { f <- function() print("hello two!") ; d <- 1:100 },
    foo2 = { f <- function() print("hello two again!") ; d <- 11:200 }
  )

  expect_output(foo1::f(), "hello two!", fixed = TRUE)
  expect_output(foo2::f(), "hello two again!", fixed = TRUE)
  expect_equal(foo1::d, 1:100)
  expect_equal(foo2::d, 11:200)

})
