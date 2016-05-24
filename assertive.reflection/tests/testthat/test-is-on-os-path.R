test_that("test.is_on_os_path.made_up_paths.returns_false_for_all", {
  paths <- c("a made up path", "path with bad chars !@#$%^&*(){}[]<>;:/?'")
  expect_false(any(is_on_os_path(paths)))
})

test_that("test.is_on_os_path.os_paths.returns_true_for_all", {
  sep <- if(.Platform$OS.type == "windows") ";" else ":"
  paths <- strsplit(Sys.getenv("PATH"), sep)[[1]]
  expect_true(all(is_on_os_path(paths)))
})

