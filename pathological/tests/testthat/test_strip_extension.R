test_that(
  "strip_extension works with paths with no directory and a single extension in the filename.",
  {
    x <- "foo.tgz"
    expected <- list(
      na    = "foo", 
      true  = file.path(getwd(), "foo", fsep = "/"), 
      false = "foo"
    )
    expected <- lapply(expected, function(y) setNames(y, x))
    expect_equal(strip_extension(x), expected$na)
    expect_equal(strip_extension(x, include_dir = TRUE), expected$true)
    expect_equal(strip_extension(x, include_dir = FALSE), expected$false)
  }
)

test_that(
  "strip_extension works with paths with a directory and a single extension in the filename.",
  {
    x <- "somedir/foo.tgz"
    expected <- list(
      na    = "somedir/foo", 
      true  = file.path(getwd(), "somedir", "foo", fsep = "/"), 
      false = "foo"
    )
    expected <- lapply(expected, function(y) setNames(y, x))
    expect_equal(strip_extension(x), expected$na)
    expect_equal(strip_extension(x, include_dir = TRUE), expected$true)
    expect_equal(strip_extension(x, include_dir = FALSE), expected$false)
  }
)

test_that(
  "strip_extension works with paths with no directory and a double extension in the filename.",
  {
    x <- "foo.tar.gz"
    expected <- list(
      na    = "foo", 
      true  = file.path(getwd(), "foo", fsep = "/"), 
      false = "foo"
    )
    expected <- lapply(expected, function(y) setNames(y, x))
    expect_equal(strip_extension(x), expected$na)
    expect_equal(strip_extension(x, include_dir = TRUE), expected$true)
    expect_equal(strip_extension(x, include_dir = FALSE), expected$false)
  }
)

test_that(
  "strip_extension works with paths with a directory and a double extension in the filename.",
  {
    x <- "somedir/foo.tar.gz"
    expected <- list(
      na    = "somedir/foo", 
      true  = file.path(getwd(), "somedir", "foo", fsep = "/"), 
      false = "foo"
    )
    expected <- lapply(expected, function(y) setNames(y, x))
    expect_equal(strip_extension(x), expected$na)
    expect_equal(strip_extension(x, include_dir = TRUE), expected$true)
    expect_equal(strip_extension(x, include_dir = FALSE), expected$false)
  }
)

test_that(
  "strip_extension works with paths with no directory and no extension in the filename.",
  {
    x <- "foo"
    expected <- list(
      na    = "foo", 
      true  = file.path(getwd(), "foo", fsep = "/"), 
      false = "foo"
    )
    expected <- lapply(expected, function(y) setNames(y, x))
    expect_equal(strip_extension(x), expected$na)
    expect_equal(strip_extension(x, include_dir = TRUE), expected$true)
    expect_equal(strip_extension(x, include_dir = FALSE), expected$false)
  }
)

test_that(
  "strip_extension works with paths with a directory and no extension in the filename.",
  {
    x <- "somedir/foo"
    expected <- list(
      na    = "somedir/foo", 
      true  = file.path(getwd(), "somedir", "foo", fsep = "/"), 
      false = "foo"
    )
    expected <- lapply(expected, function(y) setNames(y, x))
    expect_equal(strip_extension(x), expected$na)
    expect_equal(strip_extension(x, include_dir = TRUE), expected$true)
    expect_equal(strip_extension(x, include_dir = FALSE), expected$false)
  }
)

test_that(
  "strip_extension handles filenames containing a '.' and an extension.",
  {
    x <- "foo. bar.zip"
    expected <- list(
      na    = "foo. bar", 
      true  = file.path(getwd(), "foo. bar", fsep = "/"), 
      false = "foo. bar"
    )
    expected <- lapply(expected, function(y) setNames(y, x))
    expect_equal(strip_extension(x), expected$na)
    expect_equal(strip_extension(x, include_dir = TRUE), expected$true)
    expect_equal(strip_extension(x, include_dir = FALSE), expected$false)
  }
)

test_that(
  "strip_extension handles directories.",
  {
    x <- R.home()
    expected <- list(
      na    = x, 
      true  = r_home(), 
      false = ""
    )
    expected <- lapply(expected, function(y) setNames(y, x))
    expect_equal(strip_extension(x), expected$na)
    expect_equal(strip_extension(x, include_dir = TRUE), expected$true)
    expect_equal(strip_extension(x, include_dir = FALSE), expected$false)
  }
)
