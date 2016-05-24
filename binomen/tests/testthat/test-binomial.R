context("binomial")

test_that("binomial basic functionality works", {
  aa <- binomial("Poa")
  bb <- binomial("Poa", "annua")
  cc <- binomial("Poa", "annua", authority = "L.")

  expect_is(aa, "binomial")
  expect_is(bb, "binomial")
  expect_is(cc, "binomial")

  expect_named(aa, "genus")
  expect_named(bb, c("genus", "epithet"))
  expect_named(cc, c("genus", "epithet", "authority"))

  expect_equal(length(aa), 1)
  expect_equal(length(bb), 2)
  expect_equal(length(cc), 3)

  expect_equal(aa$genus, "Poa")
})

test_that("binomial fails well", {
  expect_equal(length(binomial()), 0)
  expect_error(binomial(stuff = 5), "unused argument")
})
