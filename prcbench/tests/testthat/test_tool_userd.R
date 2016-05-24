context("Tool: User-defined tool")
# Test create_example_func
#      create_usrtool
#

test_that("create_example_func", {

  func <- create_example_func()
  expect_true(is.function(func))
  expect_equal(names(formals(func)), "single_testset")

})

test_that("create_usrtool - tool_name", {
  func <- create_example_func()

  expect_error(create_usrtool("xyz", func), NA)

  expect_error(create_usrtool(1, func), "tool_name is not a string")
  expect_error(create_usrtool(c("xyz", "A"), func), "tool_name is not a string")
})

test_that("create_usrtool - func", {
  func <- create_example_func()

  expect_error(create_usrtool("xyz", func), NA)

  expect_error(create_usrtool("xyz", "A"), "func is not a function")
  expect_error(create_usrtool("xyz", function() NULL), "Invalid func")
  expect_error(create_usrtool("xyz", function(x) stop("s")), "Invalid func")
  expect_error(create_usrtool("xyz", function(x) list(1, 2, 3)), "Invalid func")
  expect_error(create_usrtool("xyz", function(x) list(x = 1,  y = 2, 3)),
               "Invalid func")
})

test_that("create_usrtool - R6", {
  func <- create_example_func()
  tool1 <- create_usrtool("xyz", func)[[1]]

  expect_true(is(tool1, "Toolxyz"))
  expect_true(is(tool1, "ToolIFBase"))
  expect_true(is(tool1, "R6"))
})

test_that("create_usrtool: calc_auc", {
  func <- create_example_func()

  expect_error(create_usrtool("xyz", func, calc_auc = TRUE), NA)
  expect_error(create_usrtool("xyz", func, calc_auc = FALSE), NA)

  expect_error(create_usrtool("xyz", func, calc_auc = 1),
               "calc_auc is not a flag")
  expect_error(create_usrtool("xyz", func, calc_auc = "TRUE"),
               "calc_auc is not a flag")

  tool1 <- create_usrtool("xyz", func)[[1]]
  expect_equal(environment(tool1$clone)$private$def_calc_auc, TRUE)

  tool2 <- create_usrtool("xyz", func, calc_auc = FALSE)[[1]]
  expect_equal(environment(tool2$clone)$private$def_calc_auc, FALSE)
})

test_that("create_usrtool: store_res", {
  func <- create_example_func()

  expect_error(create_usrtool("xyz", func, store_res = TRUE), NA)
  expect_error(create_usrtool("xyz", func, store_res = FALSE), NA)

  expect_error(create_usrtool("xyz", func, store_res = 1),
               "store_res is not a flag")
  expect_error(create_usrtool("xyz", func, store_res = "TRUE"),
               "store_res is not a flag")

  tool1 <- create_usrtool("xyz", func)[[1]]
  expect_equal(environment(tool1$clone)$private$def_store_res, TRUE)

  tool2 <- create_usrtool("xyz", func, store_res = FALSE)[[1]]
  expect_equal(environment(tool2$clone)$private$def_store_res, FALSE)
})
