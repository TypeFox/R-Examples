# TODO: test for is_solaris

test_that("test.is_bsd.any_os.returns_true_if_os_is_osx", {
  expected <- grepl("BSD", unname(Sys.info()["sysname"]))
  actual <- is_bsd()
  expect_equal(strip_attributes(actual), expected)
  if (!actual) {
    expect_equal(cause(actual), cause(not_this_os("BSD-based")))
  }
})

test_that("test.is_linux.any_os.returns_true_if_os_is_linux", {
  expected <- unname(Sys.info()["sysname"] == "Linux")
  actual <- is_linux()
  expect_equal(strip_attributes(actual), expected)
  if (!actual) {
    expect_equal(cause(actual), cause(not_this_os("Linux")))
  }
})

test_that("test.is_mac.any_os.returns_true_if_os_is_osx", {
  expected <- unname(Sys.info()["sysname"] == "Darwin")
  actual <- is_mac()
  expect_equal(strip_attributes(actual), expected)
  if (!actual) {
    expect_equal(cause(actual), cause(not_this_os("OS X")))
  }
})

test_that("test.is_solaris.any_os.returns_true_if_os_is_solaris", {
  expected <- unname(Sys.info()["sysname"] == "SunOS")
  actual <- is_solaris()
  expect_equal(strip_attributes(actual), expected)
  if (!actual) {
    expect_equal(cause(actual), cause(not_this_os("Solaris")))
  }
})

test_that("test.is_unix.any_os.returns_true_if_os_is_unix_based", {
  expected <- .Platform$OS.type == "unix"
  actual <- is_unix()
  expect_equal(strip_attributes(actual), expected)
  if (!actual) {
    expect_equal(cause(actual), cause(not_this_os("Unix-based")))
  }
})

test_that("test.is_windows.any_os.returns_true_if_os_is_windows", {
  expected <- .Platform$OS.type == "windows"
  actual <- is_windows()
  expect_equal(strip_attributes(actual), expected)
  if (!actual) {
    expect_equal(cause(actual), cause(not_this_os("Windows")))
  }
}) 

test_that(
  "test.not_this_os_returns_os_message",
  {
    os <-"foo"
    expected_cause <- noquote(
      paste0(
        "The operating system is not ",
        os,
        ". R reports it as: Sys.info()['sysname'] = ",
        Sys.info()['sysname'],
        ", .Platform$OS = ", 
        .Platform$OS,
        "."
      )
    )
    actual <- not_this_os(os)
    expect_equal(strip_attributes(actual), FALSE)
    expect_equal(cause(actual), expected_cause)
  }
)
