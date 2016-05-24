context("gisd functions")

test_that("gisd works", {
  skip_on_cran()

  sp <- c("Carpobrotus edulis", "Rosmarinus officinalis")
  aa <- gisd(sp, verbose = FALSE)

  expect_is(aa, "list")
  expect_named(aa, sp)
  expect_is(aa[[1]], "list")
  expect_is(aa[[2]], "list")
})

test_that("fails well - species not found when searching GBIF", {
  skip_on_cran()

  sp <- "asdfadsf"
  aa <- gisd(sp, verbose = FALSE)

  expect_is(aa, "list")
  expect_named(aa, sp)
  expect_equal(aa[[1]]$status, "Not in GISD")
})
