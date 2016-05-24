context("Travis Setup")

on_travis <- identical(Sys.getenv("TRAVIS"), "true")

test_that("ms is available on travis", {
  if (!on_travis) skip("Not on Travis-CI")
  expect_true(has_ms())
})


test_that("msms is available on travis", {
  if (!on_travis) skip("Not on Travis-CI")
  expect_true(has_msms())
})


test_that("seqgen is available on travis", {
  if (!on_travis) skip("Not on Travis-CI")
  expect_true(has_seqgen())
})


test_that("OmegaPlus is available on travis", {
  if (!on_travis) skip("Not on Travis-CI")
  expect_true(has_omega())
})


test_that("release questions are present", {
  release_questions()
})
