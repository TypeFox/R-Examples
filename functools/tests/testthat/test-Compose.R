library(functools)
context("Compose()")

test_that("Produces the correct output.", {
  expect_equal(Compose(Identity, Identity)(Identity), Identity(Identity))
  expect_equal(Compose("Identity", "Identity")(Identity), Identity(Identity))
  expect_equal(`%O%`(Identity, Identity)(Identity), Identity(Identity))
  expect_equal(Compose(Identity, Identity)(Identity), `%O%`(Identity, Identity)(Identity))

})

test_that("Produces the correct output type.", {
  expect_is(Compose(Identity, Identity), "function")
  expect_is(Compose("Identity", "Identity"), "function")
  expect_is(`%O%`(Identity, Identity)(Identity), "function")
})

test_that("Produces the correct errors.", {
  expect_error(Compose("a", "b"))
  expect_error(`%O%`("a", "b"))
})
