test_that(
  "replace_extension handles filenames with a single extension.",
  {
    x <- "somedir/foo.tgz"
    new_extension <- "NEW"
    expected <- list(
      na    = "somedir/foo.NEW", 
      true  = file.path(getwd(), "somedir", "foo.NEW", fsep = "/"), 
      false = "foo.NEW"
    )
    expected <- lapply(expected, function(y) setNames(y, x))
    expect_equal(
      replace_extension(x, new_extension), 
      expected$na
    )
    expect_equal(
      replace_extension(x, new_extension, include_dir = TRUE), 
      expected$true
    )
    expect_equal(
      replace_extension(x, new_extension, include_dir = FALSE), 
      expected$false
    )
  }
)

test_that(
  "replace_extension handles filenames with a double extension.",
  {
    x <- "somedir/foo.tar.gz"
    new_extension <- "NEW"
    expected <- list(
      na    = "somedir/foo.NEW", 
      true  = file.path(getwd(), "somedir", "foo.NEW", fsep = "/"), 
      false = "foo.NEW"
    )
    expected <- lapply(expected, function(y) setNames(y, x))
    expect_equal(
      replace_extension(x, new_extension), 
      expected$na
    )
    expect_equal(
      replace_extension(x, new_extension, include_dir = TRUE), 
      expected$true
    )
    expect_equal(
      replace_extension(x, new_extension, include_dir = FALSE), 
      expected$false
    )
  }
)

test_that(
  "replace_extension handles filenames with no extension.",
  {
    x <- "somedir/foo"
    new_extension <- "NEW"
    expected <- list(
      na    = "somedir/foo.NEW", 
      true  = file.path(getwd(), "somedir", "foo.NEW", fsep = "/"), 
      false = "foo.NEW"
    )
    expected <- lapply(expected, function(y) setNames(y, x))
    expect_equal(
      replace_extension(x, new_extension), 
      expected$na
    )
    expect_equal(
      replace_extension(x, new_extension, include_dir = TRUE), 
      expected$true
    )
    expect_equal(
      replace_extension(x, new_extension, include_dir = FALSE), 
      expected$false
    )
  }
)

test_that(
  "replace_extension handles directories.",
  {
    # This has to be a real directory since it is not possible to tell if
    # a non-existent 'foo' refers to a directory or filename.
    x <- R.home()
    new_extension <- "NEW"
    expected <- list(
      na    = x, 
      true  = r_home(), 
      false = ""
    )
    expected <- lapply(expected, function(y) setNames(y, x))
    error_rx <- "The directories .* have no file extensions to replace."
    actual <- list()
    expect_warning(
      actual$na <- replace_extension(x, new_extension), 
      error_rx
    )
    expect_warning(
      actual$true <- replace_extension(x, new_extension, include_dir = TRUE), 
      error_rx
    )
    expect_warning(
      actual$false <- replace_extension(x, new_extension, include_dir = FALSE), 
      error_rx
    )
    expect_equal(actual$na, expected$na)
    expect_equal(actual$true, expected$true)
    expect_equal(actual$false, expected$false)
  }
)

test_that(
  "replace_extension handles empty replacement extensions.",
  {
    x <- "somedir/foo.tgz"
    new_extension <- ""
    expected <- list(
      na    = "somedir/foo.", 
      true  = file.path(getwd(), "somedir", "foo.", fsep = "/"), 
      false = "foo."
    )
    expected <- lapply(expected, function(y) setNames(y, x))
    error_rx <- "'new_extension' is empty.  Did you want strip_extension instead?"
    actual <- list()
    expect_warning(
      actual$na <- replace_extension(x, new_extension), 
      error_rx
    )
    expect_warning(
      actual$true <- replace_extension(x, new_extension, include_dir = TRUE), 
      error_rx
    )
    expect_warning(
      actual$false <- replace_extension(x, new_extension, include_dir = FALSE), 
      error_rx
    )
    expect_equal(actual$na, expected$na)
    expect_equal(actual$true, expected$true)
    expect_equal(actual$false, expected$false)
  }
)
