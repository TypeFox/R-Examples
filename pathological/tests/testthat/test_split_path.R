# Testing expect_equal on lists doesn't give very revealign information upon
# failure, so expectations are split into class/length/value/name tests.
# Names especially need to be tested separately, since as of R3.1.1 there is a 
# bug in how names are printed in lists. See
# http://r.789695.n4.nabble.com/The-behaviour-of-setting-names-differs-between-lists-and-atomic-vectors-td4695790.html

test_that(
  "split_path works with a zero length input",
  {
    x <- character()
    actual <- split_path(x)
    expect_is(actual, "list")
    expect_equal(length(actual), 0L)
    expect_equal(unname(actual), list())
    expect_equal(names(actual), character())
  }
)

test_that(
  "split_path works with a NULL input",
{
  x <- NULL
  actual <- split_path(x)
  expect_is(actual, "list")
  expect_equal(length(actual), 0L)
  expect_equal(unname(actual), list())
  expect_equal(names(actual), character())
}
)

test_that(
  "split_path works with empty strings",
  {
    x <- ""
    expected1 <- character()
    actual <- split_path(x)
    expect_is(actual, "list")
    expect_equal(length(actual), 1L)
    expect_equal(unname(actual[[1]]), expected1)
    expect_equal(names(actual), x)
  }
)

test_that(
  "split_path works with missing values",
  {
    x <- NA
    expected1 <- NA_character_
    expect_warning(
      actual <- split_path(x),
      "Coercing .+ to class .{1,3}character.{1,3}\\."
    )
    expect_is(actual, "list")
    expect_equal(length(actual), 1L)
    expect_equal(unname(actual[[1]]), expected1)
    expect_equal(names(actual), NA_character_)
  }
)

test_that(
  "split_path works with absolute Windows paths with forward slashes.",
  {
    x <- "c:/foo/bar"
    expected1 <- c("c:", "foo", "bar")
    actual <- split_path(x)
    expect_is(actual, "list")
    expect_equal(length(actual), 1L)
    expect_equal(unname(actual[[1]]), expected1)
    expect_equal(names(actual), x)
  }
)

test_that(
  "split_path works with absolute Windows paths with back slashes.",
  {
    x <- "c:\\foo\\bar"
    expected1 <- c("c:", "foo", "bar")
    actual <- split_path(x)
    expect_is(actual, "list")
    expect_equal(length(actual), 1L)
    expect_equal(unname(actual[[1]]), expected1)
    expect_equal(names(actual), x)
  }
)

test_that(
  "split_path works with absolute Windows paths with mixed forward and back slashes.",
  {
    x <- "c:/foo\\bar"
    expected1 <- c("c:", "foo", "bar")
    actual <- split_path(x)
    expect_is(actual, "list")
    expect_equal(length(actual), 1L)
    expect_equal(unname(actual[[1]]), expected1)
    expect_equal(names(actual), x)
  }
)

test_that(
  "split_path works with absolute UNC paths with forward slashes.",
  {
    x <- "//foo/bar"
    expected1 <- c("//foo", "bar")
    actual <- split_path(x)
    expect_is(actual, "list")
    expect_equal(length(actual), 1L)
    expect_equal(unname(actual[[1]]), expected1)
    expect_equal(names(actual), x)
  }
)

test_that(
  "split_path works with absolute UNC paths with back slashes.",
  {
    x <- "\\\\foo\\bar"
    expected1 <- c("//foo", "bar")
    actual <- split_path(x)
    expect_is(actual, "list")
    expect_equal(length(actual), 1L)
    expect_equal(unname(actual[[1]]), expected1)
    expect_equal(names(actual), x)
  }
)


