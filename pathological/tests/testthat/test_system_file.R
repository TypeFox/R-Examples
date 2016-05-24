test_that(
  "system_file works correctly with no inputs",
  {
    actual <- system_file()
    expected <- normalizePath(find.package("base"), "/")
    expect_equal(actual, expected)
  }
)

test_that(
  "system_file works correctly with only a package input",
  {
    actual <- system_file(package = "stats")
    expected <- normalizePath(find.package("base"), "/")
  }
)

test_that(
  "system_file works correctly with a file inside the base package",
  {
    actual <- system_file("help", "AnIndex")
    expected <- normalizePath(
      file.path(find.package("base"), "help", "AnIndex"), 
      "/"
    )
  }
)

test_that(
  "system_file works correctly with a file inside another package",
  {
    actual <- system_file("R", "graphics.rdb", package = "graphics")
    expected <- normalizePath(
      file.path(find.package("graphics"), "R", "graphics.rdb"), 
      "/"
    )
  }
)
