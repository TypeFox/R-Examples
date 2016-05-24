context("New Argument Check Environment")

test_that("New Argument Check Environment Has Correct Classes",
{
  expect_equal(class(newArgCheck()),
               c("ArgCheck", "environment"))
})

test_that("New Argument Check Has Correct Initial Values",
{
  Check <- newArgCheck()
  expect_equal(mget(c("n_error", "error_msg", "n_warn", "warn_msg", 
                      "n_message", "message_msg"),
                    envir = Check),
               list(n_error = 0, error_msg = NULL,
                    n_warn = 0, warn_msg = NULL,
                    n_message = 0, message_msg = NULL))
})