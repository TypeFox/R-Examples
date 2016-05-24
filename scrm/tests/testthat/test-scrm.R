context("scrm")


test_that("parsing arguments works", {
  expect_error(scrm("10 -t 5"))
  expect_error(scrm("10 1 -t"))

  scrm("5 2 -t 5")
  scrm("5 1 -t 5 -r 10 100")
  scrm("5 1 -T")
  scrm("5 1 -I 2 2 3 1.5 -r 5 100 -L")

  expect_error(scrm("10 1 -oSFS"))
  scrm("10 1 -oSFS -t 2")
})


test_that("runs a reproducible", {
  set.seed(119)
  res1 <- scrm("10 1 -r 1 100 -t 2")
  set.seed(119)
  res2 <- scrm("10 1 -r 1 100 -t 2")
  expect_equal(res1, res2)
})


test_that("warning is given when using -seed", {
  expect_warning(scrm("10 1 -seed 17 -L"))
  expect_warning(scrm("10 1 -seed 1 2 3 -L"))
})


test_that("writing into a file works", {
  file <- tempfile("scrm_")
  scrm("10 1 -r 10 100 -t 5", file = file)
  expect_true(file.exists(file))
  output <- scan(file, what="character", sep="\n", quiet = TRUE)
  expect_equal(output[1], "scrm 10 1 -r 10 100 -t 5")
  expect_equal(output[2], "//")
  expect_true(output[3] != "")
  unlink(file)

  # Test relative path
  work_dir <- getwd()
  setwd(tempdir())
  file <- "scrm_test_file"
  scrm("10 1 -t 5", file = file)
  expect_true(file.exists(file))
  expect_true(file.exists(paste(tempdir(), file, sep="/")))
  unlink(file)

  # Test with "//" instead of "/"
  file <- "scrm_test_dir//scrm_test_file"
  expect_error(scrm("10 1 -t 5", file = file))
  dir.create("scrm_test_dir")
  scrm("10 1 -t 5", file = file)
  expect_true(file.exists(file))
  expect_true(file.exists(paste(tempdir(), file, sep="/")))

  # Test with "\" instead of "/"
  file <- "scrm_test_dir\\scrm_test_file2"
  scrm("10 1 -t 5", file = file)
  expect_true(file.exists(file))
  expect_true(file.exists(paste(tempdir(), file, sep="/")))
  unlink("scrm_test_dir", recursive = TRUE)
  setwd(work_dir)
})


test_that("-m & -em are available", {
  sum_stats <- scrm("10 1 -I 2 6 4 -em 1.0 2 1 0.5 -em 1.0 1 2 0.7 -T")
  expect_equal(length(sum_stats$trees), 1)
  sum_stats <- scrm("10 1 -I 2 6 4 -m 2 1 0.5 -m 1 2 0.7 -T")
  expect_equal(length(sum_stats$trees), 1)
})


test_that("printing help & verison information works", {
  for (version in list(scrm("-v"), scrm("--version"))) {
    expect_that(version, is_a("list"))
    expect_that(version$version, is_a("character"))
  }

  expect_error(scrm("-h"))
  expect_error(scrm("--help"))
})
