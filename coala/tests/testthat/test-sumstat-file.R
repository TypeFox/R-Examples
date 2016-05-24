context("SumStat File")

test_that("File statistic requires files", {
  model <- model_theta_tau()
  expect_false(requires_files(model))
  tmp_folder <- tempfile("requires_files_test")
  expect_true(requires_files(model + sumstat_file(tmp_folder)))
  unlink(tmp_folder, recursive = TRUE)
})


test_that("File statistic works", {
  folder <- tempfile("sumstat_file_test")
  stat <- sumstat_file(folder)

  files <- c(tempfile("test1"), tempfile("test2"))
  cat("test1", file = files[1])
  cat("test2", file = files[2])

  expect_equal(stat$get_name(), "file")
  files_copy <- stat$calculate(NULL, NULL, files, NULL)

  expect_true(file.exists(folder))
  expect_true(all(file.exists(files_copy)))
  expect_equal(scan(files_copy[1], what = "character", quiet = TRUE), "test1")
  expect_equal(scan(files_copy[2], what = "character", quiet = TRUE), "test2")

  unlink(c(files, folder), recursive = TRUE)
})
