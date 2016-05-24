context("argCheck Incrementation")

test_that("argCheck Incrementation",
{
  Check <- newArgCheck()
  addError("New Error",
           Check)
  addWarning("New Warning",
             Check)
  expect_equal(mget(c("n_error", "error_msg", "n_warn", "warn_msg"),
                    envir = Check),
               list(n_error = 1, error_msg = "New Error",
                    n_warn = 1, warn_msg = "New Warning"))
})