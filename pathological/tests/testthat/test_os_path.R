test_that(
  "os_path works when the OS PATH environment variable has not been set.",
  {
    old_path <- Sys.getenv("PATH")
    on.exit(Sys.setenv(PATH = old_path))
    Sys.unsetenv("PATH")
    expected <- character()
    warn <- "The 'PATH' environment variable is unset or empty."
    expect_warning(actual_std_is_false <- os_path(standardize = FALSE), warn)
    expect_warning(actual_std_is_true <- os_path(), warn)
    expect_equal(actual_std_is_false, expected)
    expect_equal(actual_std_is_true, setNames(expected, character()))
  }
)

test_that(
  "os_path works when the OS PATH environment variable is an empty string.",
  {
    old_path <- Sys.getenv("PATH")
    on.exit(Sys.setenv(PATH = old_path))
    Sys.setenv(PATH = "")
    expected <- character()
    warn <- "The 'PATH' environment variable is unset or empty."
    expect_warning(actual_std_is_false <- os_path(standardize = FALSE), warn)
    expect_warning(actual_std_is_true <- os_path(), warn)
    expect_equal(actual_std_is_false, expected)
    expect_equal(actual_std_is_true, setNames(expected, character()))
  }
)



test_that(
  "os_path works when the OS PATH environment variable is a non-empty string.",
  {
    splitter <- if(assertive.reflection::is_windows()) ";" else ":"
    if(Sys.getenv("PATH") == "")
    {
      # PATH is empty.  Need to make one up.
      Sys.setenv(PATH = paste(R.home("bin"), getwd(), collapse = splitter))
      on.exit(Sys.setenv(PATH = ""))
    }
    path <- Sys.getenv("PATH")
    expected_std_is_false <- strsplit(path, splitter)[[1]]
    expected_std_is_true <- standardize_path(expected_std_is_false, include_names = FALSE)
    expect_equal(os_path(standardize = FALSE), expected_std_is_false)
    expect_equal(os_path(), expected_std_is_true)
  }
)
