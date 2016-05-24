context("Tool: PerfMeas")
# Test ToolPerfMeas
#      create_toolset
#

test_that("ToolPerfMeas - R6ClassGenerator", {
  expect_true(is(ToolPerfMeas, "R6ClassGenerator"))
  expect_equal(attr(ToolPerfMeas, "name"), "ToolPerfMeas_generator")

  expect_equal(grep("PerfMeas",
                    body(ToolPerfMeas$private_methods$f_wrapper))[[1]], 2)
})

test_that("ToolPerfMeas - R6", {
  toolset <- ToolPerfMeas$new()

  expect_true(is(toolset, "ToolPerfMeas"))
  expect_true(is(toolset, "ToolIFBase"))
  expect_true(is(toolset, "R6"))
})

test_that("create_toolset", {
  toolset1 <- create_toolset("PERF")[[1]]
  expect_true(is(toolset1, "ToolPerfMeas"))
  expect_equal(toolset1$get_toolname(), "PerfMeas")

  toolset2 <- create_toolset("perf")[[1]]
  expect_true(is(toolset2, "ToolPerfMeas"))
  expect_equal(toolset2$get_toolname(), "PerfMeas")
})

test_that("create_toolset: calc_auc", {
  toolset1 <- create_toolset("PerfMeas")[[1]]
  expect_equal(environment(toolset1$clone)$private$def_calc_auc, TRUE)

  toolset2 <- create_toolset("PerfMeas", calc_auc = FALSE)[[1]]
  expect_equal(environment(toolset2$clone)$private$def_calc_auc, FALSE)
})

test_that("create_toolset: store_res", {
  toolset1 <- create_toolset("PerfMeas")[[1]]
  expect_equal(environment(toolset1$clone)$private$def_store_res, TRUE)

  toolset2 <- create_toolset("PerfMeas", store_res = FALSE)[[1]]
  expect_equal(environment(toolset2$clone)$private$def_store_res, FALSE)
})
