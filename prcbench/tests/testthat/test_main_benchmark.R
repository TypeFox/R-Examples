context("Main: Benchmark")
# Test run_benchmark
#

.check_microbenchmark <- function() {
  use_sys_time <- FALSE
  if (!requireNamespace("microbenchmark", quietly = TRUE)) {
    use_sys_time <- TRUE
  }

  use_sys_time
}

test_that("run_benchmark", {
  tfunc <- function(use_sys_time) {
    testset <- create_testset("bench", c("b10", "i10"))
    toolset <- create_toolset(set_names = c("crv5", "auc5"))
    res1 <- run_benchmark(testset, toolset, times = 2,
                          use_sys_time = use_sys_time)[["tab"]]

    num_tools <- 5
    num_toolsets <- 2
    num_datasets <- 2
    expect_equal(nrow(res1), num_tools * num_toolsets * num_datasets)
    expect_equal(length(unique(res1$toolname)), num_tools)
    expect_equal(length(unique(res1$testset)), num_datasets)
    expect_equal(length(unique(res1$toolset)), num_toolsets)
    expect_true(all(res1$neval == 2))
  }

  use_sys_time <- .check_microbenchmark()
  tfunc(use_sys_time)

  if (!use_sys_time) {
    tfunc(TRUE)
  }
})

test_that("run_benchmark: testset", {
  tfunc <- function(use_sys_time) {
    testset <- create_testset("bench", c("b10", "i10"))
    toolset <- create_toolset(set_names = c("crv5", "auc5"))

    expect_error(run_benchmark(testset, toolset, use_sys_time = use_sys_time),
                 NA)

    expect_error(run_benchmark(1, toolset, use_sys_time = use_sys_time),
                 "testset is not a list")
    expect_error(run_benchmark("1", toolset, use_sys_time = use_sys_time),
                 "testset is not a list")
    expect_error(run_benchmark(list(), toolset, use_sys_time = use_sys_time),
                 "not greater than 0")
    expect_error(run_benchmark(toolset, toolset, use_sys_time = use_sys_time),
                 "Invalid testset")
  }

  use_sys_time <- .check_microbenchmark()
  tfunc(use_sys_time)

  if (!use_sys_time) {
    tfunc(TRUE)
  }
})

test_that("run_benchmark: toolset", {
  tfunc <- function(use_sys_time) {
    testset <- create_testset("bench", c("b10", "i10"))
    toolset <- create_toolset(set_names = c("crv5", "auc5"))

    expect_error(run_benchmark(testset, toolset, use_sys_time = use_sys_time),
                 NA)

    expect_error(run_benchmark(testset, 1, use_sys_time = use_sys_time),
                 "toolset is not a list")
    expect_error(run_benchmark(testset, "1", use_sys_time = use_sys_time),
                 "toolset is not a list")
    expect_error(run_benchmark(testset, list(), use_sys_time = use_sys_time),
                 "not greater than 0")
    expect_error(run_benchmark(testset, testset, use_sys_time = use_sys_time),
                 "Invalid toolset")
  }

  use_sys_time <- .check_microbenchmark()
  tfunc(use_sys_time)

  if (!use_sys_time) {
    tfunc(TRUE)
  }
})

test_that("run_benchmark: times", {
  tfunc <- function(use_sys_time) {
    testset <- create_testset("bench", "b10")
    toolset <- create_toolset(set_names = "def5")

    res <- run_benchmark(testset, toolset, times = 1,
                         use_sys_time = use_sys_time)[["tab"]]
    expect_true(all(res$neval == 1))

    expect_error(run_benchmark(testset, toolset, times = 0,
                               use_sys_time = use_sys_time),
                 "times not greater than 0")
    expect_error(run_benchmark(testset, toolset, times = "1",
                               use_sys_time = use_sys_time),
                 "times is not a number")
  }

  use_sys_time <- .check_microbenchmark()
  tfunc(use_sys_time)

  if (!use_sys_time) {
    tfunc(TRUE)
  }
})

test_that("run_benchmark: unit", {
  tfunc <- function(use_sys_time) {
    testset <- create_testset("bench", "b10")
    toolset <- create_toolset(set_names = "def5")

    expect_error(run_benchmark(testset, toolset, unit = "ns",
                               use_sys_time = use_sys_time), NA)
    expect_error(run_benchmark(testset, toolset, unit = "us",
                               use_sys_time = use_sys_time), NA)
    expect_error(run_benchmark(testset, toolset, unit = "ms",
                               use_sys_time = use_sys_time), NA)
    expect_error(run_benchmark(testset, toolset, unit = "s",
                               use_sys_time = use_sys_time), NA)

    expect_error(run_benchmark(testset, toolset, unit = "eps",
                               use_sys_time = use_sys_time), NA)
    expect_error(run_benchmark(testset, toolset, unit = "relative",
                               use_sys_time = use_sys_time), NA)

    expect_error(run_benchmark(testset, toolset, unit = "ss",
                               use_sys_time = use_sys_time),
                 "is not TRUE")
    expect_error(run_benchmark(testset, toolset, unit = 1,
                               use_sys_time = use_sys_time),
                 "unit is not a string")
  }

  use_sys_time <- .check_microbenchmark()
  tfunc(use_sys_time)

  if (!use_sys_time) {
    tfunc(TRUE)
  }
})