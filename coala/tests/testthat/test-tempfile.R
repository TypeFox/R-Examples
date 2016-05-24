context("Tempfiles")

test_that("Threadsafe tempfile implementation works", {
  tmp <- tempfile("blub")
  expect_equal(grep("coala", tmp), 1)
  expect_equal(grep(Sys.getpid(), tmp), 1)
  expect_equal(grep("blub", tmp), 1)
})


test_that("there are no remaining temporary files", {
  temp_files <- list.files(tempdir(), pattern = "^coala-[0-9]+-")
  for (temp_file in temp_files) {
    cat("temp file not deleted: ", temp_file, "\n")
    warning("temp file not deleted: ", temp_file)
  }
  expect_equal(length(temp_files), 0)
})


test_that("no seedms files were created", {
  expect_false(file.exists("seedms"))
})
