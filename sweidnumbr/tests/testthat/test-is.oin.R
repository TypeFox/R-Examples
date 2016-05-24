
test_oin <- c("556000-4615", "232100-0156", "802002-4280", "8020024280", "802002A4280", "8020A2-4280","801002-4280", "801002-428 ")

context("is.oin")

test_that(desc="is.pin",{
  expect_is(is.oin(oin = test_oin), "logical")
  expect_equal(is.oin(oin = test_oin), expected = FALSE)
  expect_equal(is.oin(oin = suppressWarnings(as.oin(test_oin))), expected = TRUE)
})
