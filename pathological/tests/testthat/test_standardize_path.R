test_that(
  "standardize_path works with a zero length input",
  {
    x <- character()
    x2 <- NULL
    expected <- setNames(character(), character())
    expect_equal(standardize_path(x), expected)
    expect_equal(standardize_path(x2), expected)
  }
)

test_that(
  "standardize_path works with empty strings",
  {
    x <- ""
    expected <- setNames("", "")
    expect_equal(standardize_path(x), expected)
  }
)

test_that(
  "standardize_path works with missing inputs",
  {
    x <- NA
    expected <- setNames(NA_character_, NA)
    expect_warning(
      actual <- standardize_path(x),
      "Coercing .+ to class .{1,3}character.{1,3}\\."
    )
    expect_equal(actual, expected)
  }
)

test_that(
  "standardize_path works with relative paths with forward slashes.",
  {
    x <- "somedir/foo.tgz"
    pwd <- getwd()
    expected <- setNames(
      file.path(pwd, "somedir", "foo.tgz", fsep = "/"),
      x
    )
    expect_equal(standardize_path(x), expected)
  }
)

test_that(
  "standardize_path works with relative paths with back slashes.",
  {
    x <- "somedir\\foo.tgz"
    pwd <- getwd()
    expected <- setNames(
      file.path(pwd, "somedir", "foo.tgz", fsep = "/"),
      x
    )
    expect_equal(standardize_path(x), expected)
  }
)

test_that(
  "standardize_path works with relative paths with mixed forward and back slashes.",
  {
    x <- "somedir/another dir\\foo.tgz"
    pwd <- getwd()
    expected <- setNames(
      file.path(pwd, "somedir", "another dir", "foo.tgz", fsep = "/"),
      x
    )
    expect_equal(standardize_path(x), expected)
  }
)

test_that(
  "standardize_path works with absolute Windows paths with forward slashes.",
  {
    x <- "c:/foo/bar"
    expected <- setNames("c:/foo/bar", x)
    expect_equal(standardize_path(x), expected)
  }
)

test_that(
  "standardize_path works with absolute Windows paths with back slashes.",
  {
    x <- "c:\\foo\\bar"
    expected <- setNames("c:/foo/bar", x)
    expect_equal(standardize_path(x), expected)
  }
)

test_that(
  "standardize_path works with absolute Windows paths with mixed forward and back slashes.",
  {
    x <- "c:/foo\\bar"
    expected <- setNames("c:/foo/bar", x)
    expect_equal(standardize_path(x), expected)
  }
)

test_that(
  "standardize_path works with absolute UNC paths with forward slashes.",
  {
    x <- "//foo/bar"
    expected <- setNames("//foo/bar", x)
    expect_equal(standardize_path(x), expected)
  }
)

test_that(
  "standardize_path works with absolute UNC paths with back slashes.",
  {
    x <- "\\\\foo\\bar"
    expected <- setNames("//foo/bar", x)
    expect_equal(standardize_path(x), expected)
  }
)
