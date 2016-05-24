context("requirePackages")

test_that("requirePackages", {
  expect_equal(requirePackages("base"), c(base=TRUE))
  expect_equal(requirePackages("xxx", stop = FALSE, suppress.warnings = TRUE), c(xxx=FALSE))
  expect_error(requirePackages("xxx", suppress.warnings=TRUE), "Please install the following packages: xxx")
  expect_equal(requirePackages(c("xxx", "base"), stop=FALSE, suppress.warnings=TRUE), c(xxx=FALSE, base=TRUE))
  expect_equal(requirePackages(c("base", "xxx"), stop=FALSE, suppress.warnings=TRUE), c(base=TRUE, xxx=FALSE))
  expect_error(requirePackages(c("base", "xxx"), suppress.warnings=TRUE), "Please install the following packages: xxx")
  expect_error(requirePackages(c("base", "xxx"), why="test", suppress.warnings=TRUE), "For test please install the following packages: xxx")

  # test loading vs. attaching using the codetools package
  expect_equal(requirePackages("codetools", default.method = "load"), c(codetools=TRUE))
  expect_true("codetools" %in% loadedNamespaces())
  expect_false("package:codetools" %in% search())
  expect_equal(requirePackages("!codetools", default.method = "load"), c(codetools=TRUE))
  expect_true("package:codetools" %in% search())
})
