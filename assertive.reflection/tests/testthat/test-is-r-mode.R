test_that("test.is_batch_mode.any_mode.returns_true_if_called_from_batch_mode", 
  {
    expected <- !is.na(Sys.getenv("R_BATCH", NA))
    actual <- is_batch_mode()
    expect_equal(strip_attributes(actual), expected)
    if (!actual) {
      expect_equal(cause(actual), noquote("R is not running in batch mode."))
    }
  })

test_that("test.is_interactive.any_mode.returns_true_if_r_runs_interactively", 
{
  expected <- interactive()
  actual <- is_interactive()
  expect_equal(strip_attributes(actual), expected)
  if (!actual) {
    expect_equal(cause(actual), noquote("R is not running interactively."))
  }
})

test_that("test.is_r.r_or_s.returns_true_if_is_r", {
  expected <- exists("is.R") && is.function(is.R) && is.R()
  actual <- is_r()
  expect_equal(strip_attributes(actual), expected)
  if (!actual) {
    expect_equal(cause(actual), noquote("You are not running R."))
  }
})


test_that("test.is_r_slave.any_os.returns_true_if_slave_in_command_args", {
  expected <- "--slave" %in% commandArgs()
  actual <- is_r_slave()
  expect_equal(strip_attributes(actual), expected)
  if (!actual) {
    expect_equal(
      cause(actual), 
      noquote("You are not running a slave instance of R.")
    )
  }
})
