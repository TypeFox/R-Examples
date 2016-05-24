
context("oin_ctrl")

test_that(desc="oin_ctrl",{
  test_oin <- c("556000-4615", "232100-0156", "802002-4280", "556000-4617", "232100-0152", "802002-4281")
  test_oin_res <- rep(FALSE, 6)
  test_oin_res[1:3] <- rep(TRUE, 3)
  
  expect_is(oin_ctrl(oin = test_oin), "logical")
  expect_equal(oin_ctrl(oin = test_oin), expected = test_oin_res)
})


test_that(desc="Expect force_logical",{
  num_to_check <- c("202100-6255","121212-1212","19121212-1212","121212+1212","1212121212",
                    1212121212, NA, Inf, TRUE, F, "foo", 123, 456L)
  suppressWarnings(expect_equal(oin_ctrl(num_to_check), 
               c(TRUE, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA)))
  expect_equal(oin_ctrl(num_to_check, force_logical=TRUE), 
               c(TRUE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE))
})


test_that(desc="Expect that NA don't cause error",{
  num_to_check <- c("202100-6255","121212-1212","19121212-1212","121212+1212","1212121212",
                    1212121212, NA, Inf, TRUE, F, "foo", 123, 456L)
  expect_silent(suppressWarnings(oin_ctrl(num_to_check)))
})

test_that(desc="Expect wrong input lenght",{
  test_oin <- c("556000-4615", "556000-461", "556000-46", "556000-4", "556000-", "556000","556000-46155")
  test_oin_res <- rep(NA, 7)
  test_oin_res[1] <- TRUE
  test_oin_res_logi <- rep(FALSE, 7)
  test_oin_res_logi[1] <- TRUE

  suppressWarnings(expect_equal(oin_ctrl(oin = test_oin), expected = test_oin_res))
  expect_is(oin_ctrl(oin = test_oin, force_logical=TRUE), "logical")
  expect_equal(oin_ctrl(oin = test_oin, force_logical=TRUE), expected = test_oin_res_logi)
})
