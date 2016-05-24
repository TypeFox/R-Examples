context("Tool: PRROC")
# Test ToolPRROC
#      create_toolset
#

test_that("ToolPRROC - R6ClassGenerator", {
  expect_true(is(ToolPRROC, "R6ClassGenerator"))
  expect_equal(attr(ToolPRROC, "name"), "ToolPRROC_generator")

  expect_true(is.function(ToolPRROC$public_methods$set_curve))
  expect_true(is.function(ToolPRROC$public_methods$set_minStepSize))

  expect_equal(grep(".prroc_wrapper",
                    body(ToolPRROC$private_methods$f_wrapper))[[1]], 2)
})

test_that("ToolPRROC - R6", {
  toolset <- ToolPRROC$new()

  expect_true(is(toolset, "ToolPRROC"))
  expect_true(is(toolset, "ToolIFBase"))
  expect_true(is(toolset, "R6"))

  expect_true(is.function(toolset[["set_curve"]]))
  expect_true(is.function(toolset[["set_minStepSize"]]))
})

test_that("ToolPRROC$new(curve)", {
  toolset1 <- ToolPRROC$new()
  expect_equal(environment(toolset1$clone)$private$curve, TRUE)

  toolset2 <- ToolPRROC$new(curve = FALSE)
  expect_equal(environment(toolset2$clone)$private$curve, FALSE)
})

test_that("ToolPRROC$new(minStepSize)", {
  toolset1 <- ToolPRROC$new()
  expect_equal(environment(toolset1$clone)$private$minStepSize, 0.01)

  toolset2 <- ToolPRROC$new(minStepSize = 0.05)
  expect_equal(environment(toolset2$clone)$private$minStepSize, 0.05)
})

test_that("create_toolset", {
  toolset1 <- create_toolset("PRR")[[1]]
  expect_true(is(toolset1, "ToolPRROC"))
  expect_equal(toolset1$get_toolname(), "PRROC")

  toolset2 <- create_toolset("prr")[[1]]
  expect_true(is(toolset2, "ToolPRROC"))
  expect_equal(toolset2$get_toolname(), "PRROC")
})

test_that("create_toolset: calc_auc", {
  toolset1 <- create_toolset("PRROC")[[1]]
  expect_equal(environment(toolset1$clone)$private$def_calc_auc, TRUE)

  toolset2 <- create_toolset("PRROC", calc_auc = FALSE)[[1]]
  expect_equal(environment(toolset2$clone)$private$def_calc_auc, FALSE)
})

test_that("create_toolset: store_res", {
  toolset1 <- create_toolset("PRROC")[[1]]
  expect_equal(environment(toolset1$clone)$private$def_store_res, TRUE)

  toolset2 <- create_toolset("PRROC", store_res = FALSE)[[1]]
  expect_equal(environment(toolset2$clone)$private$def_store_res, FALSE)
})
