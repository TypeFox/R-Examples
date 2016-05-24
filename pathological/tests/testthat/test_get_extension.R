test_that(
  "get_extension works with paths with no directory and a single extension in the filename.",
  {
    x <- "foo.tgz"
    expected <- "tgz"
    names(expected) <- x
    expect_equal(get_extension(x), expected)
  }
)

test_that(
  "get_extension works with paths with a directory and a single extension in the filename.",
  {
    x <- "somedir/foo.tgz"
    expected <- "tgz"
    names(expected) <- x
    expect_equal(get_extension(x), expected)
  }
)

test_that(
  "get_extension works with paths with no directory and a double extension in the filename.",
  {
    x <- "foo.tar.gz"
    expected <- "tar.gz"
    names(expected) <- x
    expect_equal(get_extension(x), expected)
  }
)

test_that(
  "get_extension works with paths with a directory and a double extension in the filename.",
  {
    x <- "somedir/foo.tar.gz"
    expected <- "tar.gz"
    names(expected) <- x
    expect_equal(get_extension(x), expected)
  }
)

test_that(
  "get_extension works with paths with no directory and no extension in the filename.",
  {
    x <- "foo"
    expected <- ""
    names(expected) <- x
    expect_equal(get_extension(x), expected)
  }
)

test_that(
  "get_extension works with paths with a directory and no extension in the filename.",
  {
    x <- "somedir/foo"
    expected <- ""
    names(expected) <- x
    expect_equal(get_extension(x), expected)
  }
)

test_that(
  "get_extension handles filenames containing a '.' and an extension.",
  {
    x <- "foo. bar.zip"
    expected <- "zip"
    names(expected) <- x
    expect_equal(get_extension(x), expected)
  }
)

test_that(
  "get_extension handles directories.",
  {
    x <- R.home()
    expected <- ""
    names(expected) <- x
    expect_equal(get_extension(x), expected)
  }
)
