context("Tool: IFBase")
# Test ToolIFBase
#

test_that("ToolIFBase - R6ClassGenerator", {
  expect_true(is(ToolIFBase, "R6ClassGenerator"))
  expect_equal(attr(ToolIFBase, "name"), "ToolIFBase_generator")

  expect_true(is.function(ToolIFBase$public_methods$call))
  expect_true(is.function(ToolIFBase$public_methods$get_toolname))
  expect_true(is.function(ToolIFBase$public_methods$get_setname))
  expect_true(is.function(ToolIFBase$public_methods$get_result))
  expect_true(is.function(ToolIFBase$public_methods$get_x))
  expect_true(is.function(ToolIFBase$public_methods$get_y))
  expect_true(is.function(ToolIFBase$public_methods$get_auc))
  expect_true(is.function(ToolIFBase$public_methods$print))
})

test_that("ToolIFBase - R6", {
  tool_obj <- ToolIFBase$new()

  expect_true(is(tool_obj, "ToolIFBase"))
  expect_true(is(tool_obj, "R6"))

  expect_true(is.function(tool_obj[["call"]]))
  expect_true(is.function(tool_obj[["get_toolname"]]))
  expect_true(is.function(tool_obj[["get_setname"]]))
  expect_true(is.function(tool_obj[["get_result"]]))
  expect_true(is.function(tool_obj[["get_x"]]))
  expect_true(is.function(tool_obj[["get_y"]]))
  expect_true(is.function(tool_obj[["get_auc"]]))
  expect_true(is.function(tool_obj[["print"]]))
})

test_that("ToolIFBase$get_x", {
  tool_obj <- ToolIFBase$new()

  expect_true(is.na(tool_obj$get_x()))
})

test_that("ToolIFBase$get_y", {
  tool_obj <- ToolIFBase$new()

  expect_true(is.na(tool_obj$get_y()))
})

test_that("ToolIFBase$get_auc", {
  tool_obj <- ToolIFBase$new()

  expect_true(is.na(tool_obj$get_auc()))
})
