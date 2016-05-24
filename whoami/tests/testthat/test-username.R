
context("Username")

test_that("username works", {

  user <- Sys.getenv('LOGNAME')
  on.exit(Sys.setenv(LOGNAME = user), add = TRUE)
  Sys.setenv(LOGNAME = 'jambajoe')

  expect_equal(username(), 'jambajoe')

})

test_that("username fallback works", {

  LOGNAME  <- Sys.getenv("LOGNAME")
  USER     <- Sys.getenv("USER")
  LNAME    <- Sys.getenv("LNAME")
  USERNAME <- Sys.getenv("USERNAME")
  on.exit(Sys.setenv(LOGNAME = LOGNAME, USER = USER,
                     LNAME = LNAME, USERNAME = USERNAME), add = TRUE)

  expect_match(username(), ".*")
  
})
