context("more.aa")

test_that("unmatched proteins produce a message and NA row", {
  expect_message(more.aa(c("ACCA", "XXX"), "Eco"), "XXX was not matched")
  expect_equal(nrow(more.aa(c("ACCA", "XXX"), "Eco")), 2)
})

