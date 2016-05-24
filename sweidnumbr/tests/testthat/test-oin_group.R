
test_oin <- c("556000-4615", "232100-0156", "802002-4280")
test_oin_res <- c("Aktiebolag", "Stat, landsting, kommuner, f\u00F6rsamlingar", "Ideella f\u00F6reningar och stiftelser")

context("oin_group")

test_that(desc="oin_group",{
  expect_is(suppressWarnings(oin_group(oin = test_oin)), "factor")
  expect_equal(as.character(oin_group(oin = test_oin)), expected = test_oin_res)
})
