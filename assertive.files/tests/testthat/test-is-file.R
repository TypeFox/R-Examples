test_that("test.is_dir.some_paths.returns_true_when_path_is_dir", 
{
  file.create(tmp <- tempfile())
  on.exit(unlink(tmp))
  x <- c(
    R.home(), 
    tmp, 
    "~not a real directory~"
  )
  expected <- c(TRUE, FALSE, FALSE)
  expect_equal(
    assertive.base::strip_attributes(actual <- is_dir(x)), 
    expected
  )
  expect_equal(names(actual), as.character(x))
  expect_equal(
    cause(actual),
    noquote(c("", "file", "nonexistent"))
  )
})

test_that("test.is_executable_file.r_exes.returns_true", {
  x <- dir(R.home("bin"), "\\.exe$", full.names = TRUE)
  expected <- rep.int(TRUE, length(x))
  names(expected) <- x
  expect_equal(suppressWarnings(is_executable_file(x)), expected)
  if(.Platform$OS.type == "windows")
  {
    expect_warning(
      is_executable_file(x), 
      "This function depends on file.access, which can give unexpected results under Windows."
    )
  }
})

test_that("expect_named(actual)", 
{
  tf <- tempfile()
  file.create(tf)
  on.exit(unlink(tf))
  x <- c("~", getwd(), tf, "~not an existing file~", NA)
  expected <- c(TRUE, TRUE, TRUE, FALSE, FALSE)
  expect_equal(
    assertive.base::strip_attributes(actual <- is_existing_file(x)), 
    expected
  )
  expect_named(actual)
  expect_equal(
    cause(actual),
    noquote(rep.int(c("", "nonexistent"), c(3, 2)))
  )
})

test_that("test.is_library.some_paths.returns_true_when_path_is_library", {
  x <- c(.libPaths(), "a made up directory")
  n_libs <- length(.libPaths())
  expected <- rep.int(c(TRUE, FALSE), c(n_libs, 1))
  expect_equal(
    assertive.base::strip_attributes(actual <- is_library(x)), 
    expected
  )
  expect_named(actual)
  expect_equal(
    cause(actual),
    noquote(rep.int(c("", "not a lib"), c(n_libs, 1)))
  )
})

test_that("test.is_readable_file.r_bin_files.returns_true", {
  x <- dir(R.home("bin"), full.names = TRUE)
  expected <- rep.int(TRUE, length(x))
  names(expected) <- x
  expect_equal(suppressWarnings(is_readable_file(x)), expected)
  if(.Platform$OS.type == "windows")
  {
    expect_warning(
      is_readable_file(x), 
      "This function depends on file.access, which can give unexpected results under Windows."
    )
  }
})

test_that("test.is_writable_file.tempfile.returns_true", {
  file.create(x <- tempfile())
  on.exit(unlink(x))
  expected <- TRUE
  names(expected) <- x
  expect_equal(suppressWarnings(is_writable_file(x)), expected)
  if(.Platform$OS.type == "windows")
  {
    expect_warning(
      is_writable_file(x), 
      "This function depends on file.access, which can give unexpected results under Windows."
    )
  }
}) 
