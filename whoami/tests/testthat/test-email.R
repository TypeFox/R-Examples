
context("Email address")

test_that("Email address works", {

  with_mock(
    `base::system` = function(...) "jambajoe@joe.joe",
    expect_equal(email_address(), "jambajoe@joe.joe")
  )
  
})
