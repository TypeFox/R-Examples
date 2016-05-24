
context("Fallbacks")

test_that("username() falls back", {

  with_mock(
    `base::Sys.getenv` = function(...) NULL,
    `base::system` = function(...) stop(),
    u <- username(fallback = "foobar")
  )

  expect_equal(u, "foobar")
})

test_that("fullname() falls back", {

  with_mock(
    `base::system` = function(...) stop(),
    f <- fullname(fallback = "Foo Bar")
  )

  expect_equal(f, "Foo Bar")
})

test_that("email_address() falls back", {

  with_mock(
    `base::system` = function(...) stop(),
    e <- email_address(fallback = "foo@bar")
  )

  expect_equal(e, "foo@bar")
})

test_that("gh_username() falls back", {

  with_mock(
    `whoami::email_address` = function(...) "not an email",
    gh <- gh_username(fallback = "foobar")
  )

  expect_equal(gh, "foobar")

  with_mock(
    `whoami::email_address` = function(...) stop(),
    gh <- gh_username(fallback = "foobar2")
  )

  expect_equal(gh, "foobar2")
})
