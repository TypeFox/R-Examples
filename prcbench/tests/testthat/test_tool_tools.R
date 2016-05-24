context("Tool: Create wrapper objects")
# Test create_toolset
#      .rename_tool_names
#

test_that("create_toolset: tool_names", {
  expect_error(create_toolset("ROCR"), NA)
  expect_error(create_toolset("AUCCalculator"), NA)
  expect_error(create_toolset("PerfMeas"), NA)
  expect_error(create_toolset("PRROC"), NA)
  expect_error(create_toolset("precrec"), NA)
  expect_error(create_toolset(c("ROCR", "PRROC")), NA)

  expect_error(create_toolset("CROC"), "Invalid tool_names")
  expect_error(create_toolset(c("ROCR", "CROC")), "Invalid tool_names")
})

test_that("create_toolset: set_names", {
  expect_error(create_toolset(set_names = "def5"), NA)
  expect_error(create_toolset(set_names = "auc5"), NA)
  expect_error(create_toolset(set_names = "crv5"), NA)
  expect_error(create_toolset(set_names = "def4"), NA)
  expect_error(create_toolset(set_names = "auc4"), NA)
  expect_error(create_toolset(set_names = "crv4"), NA)
  expect_error(create_toolset(set_names = c("auc5", "crv4")), NA)

  expect_error(create_toolset(set_names = "crv3"), "Invalid set_names")
  expect_error(create_toolset(set_names = c("auc5", "crv3")),
                              "Invalid set_names")
})

test_that("create_toolset: calc_auc", {
  expect_error(create_toolset("ROCR", calc_auc = TRUE), NA)
  expect_error(create_toolset("ROCR", calc_auc = FALSE), NA)

  expect_error(create_toolset("ROCR", calc_auc = 1), "calc_auc is not a flag")
  expect_error(create_toolset("ROCR", calc_auc = "TRUE"),
               "calc_auc is not a flag")
})

test_that("create_toolset: store_res", {
  expect_error(create_toolset("ROCR", store_res = TRUE), NA)
  expect_error(create_toolset("ROCR", store_res = FALSE), NA)

  expect_error(create_toolset("ROCR", store_res = 1), "store_res is not a flag")
  expect_error(create_toolset("ROCR", store_res = "TRUE"),
               "store_res is not a flag")
})

test_that("create_toolset: crv5", {
  toolset1 <- create_toolset(set_names = "crv5")
  expect_equal(length(toolset1), 5)
  for (i in 1:5) {
    expect_true(is(toolset1[[i]], "R6"))
    expect_true(is(toolset1[[i]], "ToolIFBase"))
  }

  expect_true(is(toolset1[[1]], "ToolROCR"))
  expect_true(is(toolset1[[2]], "ToolAUCCalculator"))
  expect_true(is(toolset1[[3]], "ToolPerfMeas"))
  expect_true(is(toolset1[[4]], "ToolPRROC"))
  expect_true(is(toolset1[[5]], "Toolprecrec"))

  expect_true(environment(toolset1[[4]]$clone)$private$curve)

})

test_that("create_toolset: auc5", {
  toolset2 <- create_toolset(set_names = "auc5")
  expect_equal(length(toolset2), 5)
  for (i in 1:5) {
    expect_true(is(toolset2[[i]], "R6"))
    expect_true(is(toolset2[[i]], "ToolIFBase"))
  }

  expect_true(is(toolset2[[1]], "ToolROCR"))
  expect_true(is(toolset2[[2]], "ToolAUCCalculator"))
  expect_true(is(toolset2[[3]], "ToolPerfMeas"))
  expect_true(is(toolset2[[4]], "ToolPRROC"))
  expect_true(is(toolset2[[5]], "Toolprecrec"))

  expect_true(!environment(toolset2[[4]]$clone)$private$curve)
})

test_that("create_toolset: def5", {
  toolset3 <- create_toolset(set_names = "def5")
  expect_equal(length(toolset3), 5)
  for (i in 1:5) {
    expect_true(is(toolset3[[i]], "R6"))
    expect_true(is(toolset3[[i]], "ToolIFBase"))
  }

  expect_true(is(toolset3[[1]], "ToolROCR"))
  expect_true(is(toolset3[[2]], "ToolAUCCalculator"))
  expect_true(is(toolset3[[3]], "ToolPerfMeas"))
  expect_true(is(toolset3[[4]], "ToolPRROC"))
  expect_true(is(toolset3[[5]], "Toolprecrec"))
})

test_that("create_toolset: crv4", {
  toolset1 <- create_toolset(set_names = "crv4")
  expect_equal(length(toolset1), 4)
  for (i in 1:4) {
    expect_true(is(toolset1[[i]], "R6"))
    expect_true(is(toolset1[[i]], "ToolIFBase"))
  }

  expect_true(is(toolset1[[1]], "ToolROCR"))
  expect_true(is(toolset1[[2]], "ToolAUCCalculator"))
  expect_true(is(toolset1[[3]], "ToolPerfMeas"))
  expect_true(is(toolset1[[4]], "Toolprecrec"))
})

test_that("create_toolset: auc4", {
  toolset2 <- create_toolset(set_names = "auc4")
  expect_equal(length(toolset2), 4)
  for (i in 1:4) {
    expect_true(is(toolset2[[i]], "R6"))
    expect_true(is(toolset2[[i]], "ToolIFBase"))
  }

  expect_true(is(toolset2[[1]], "ToolROCR"))
  expect_true(is(toolset2[[2]], "ToolAUCCalculator"))
  expect_true(is(toolset2[[3]], "ToolPerfMeas"))
  expect_true(is(toolset2[[4]], "Toolprecrec"))
})

test_that("create_toolset: def4", {
  toolset3 <- create_toolset(set_names = "def4")
  expect_equal(length(toolset3), 4)
  for (i in 1:4) {
    expect_true(is(toolset3[[i]], "R6"))
    expect_true(is(toolset3[[i]], "ToolIFBase"))
  }

  expect_true(is(toolset3[[1]], "ToolROCR"))
  expect_true(is(toolset3[[2]], "ToolAUCCalculator"))
  expect_true(is(toolset3[[3]], "ToolPerfMeas"))
  expect_true(is(toolset3[[4]], "Toolprecrec"))
})

test_that(".rename_tool_names", {
  renamed1 <- .rename_tool_names(c("1", "2", "3", "4"))
  expect_equal(renamed1, c("1", "2", "3", "4"))

  renamed2 <- .rename_tool_names(c("1", "2", "1", "4"))
  expect_equal(renamed2, c("1", "2", "1.2", "4"))

  renamed3 <- .rename_tool_names(c("1", "2", "1", "1"))
  expect_equal(renamed3, c("1", "2", "1.2", "1.3"))

  renamed4 <- .rename_tool_names(c("1", "2", "1", "2"))
  expect_equal(renamed4, c("1", "2", "1.2", "2.2"))
})

test_that("create_toolset: single tool", {
  tool1 <- create_toolset("ROCR")
  expect_true(is(tool1[[1]], "ToolROCR"))

  tool2 <- create_toolset("AUCCalculator")
  expect_true(is(tool2[[1]], "ToolAUCCalculator"))

  tool3 <- create_toolset("PerfMeas")
  expect_true(is(tool3[[1]], "ToolPerfMeas"))

  tool4 <- create_toolset("PRROC")
  expect_true(is(tool4[[1]], "ToolPRROC"))

  tool5 <- create_toolset("precrec")
  expect_true(is(tool5[[1]], "Toolprecrec"))
})

test_that("create_toolset: multiple tools", {
  tool1 <- create_toolset(c("ROCR", "PRROC", "PerfMeas", "precrec"))
  expect_equal(names(tool1), c("ROCR", "PRROC", "PerfMeas", "precrec"))

  tool2 <- create_toolset(c("roc", "prr", "perf", "prec"))
  expect_equal(names(tool2), c("ROCR", "PRROC", "PerfMeas", "precrec"))
})

test_that("Duplicated names", {
  tool1 <- create_toolset(c("ROCR", "PRROC", "PerfMeas", "precrec"))
  expect_equal(names(tool1), c("ROCR", "PRROC", "PerfMeas", "precrec"))

  tool2 <- create_toolset(c("ROCR", "PRROC", "ROCR", "precrec"))
  expect_equal(names(tool2), c("ROCR", "PRROC", "ROCR.2", "precrec"))

  tool3 <- create_toolset(c("ROCR", "PRROC", "ROCR", "PRROC"))
  expect_equal(names(tool3), c("ROCR", "PRROC", "ROCR.2", "PRROC.2"))

  tool4 <- create_toolset(c("ROCR", "PRROC", "ROCR", "ROCR"))
  expect_equal(names(tool4), c("ROCR", "PRROC", "ROCR.2", "ROCR.3"))
})
