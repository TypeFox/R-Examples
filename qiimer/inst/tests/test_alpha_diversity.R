context('read_qiime_rarefaction')

r <- read_qiime_rarefaction(
  system.file("testdata", "observed_species.txt", package="qiimer"))

test_that("Sample IDs are parsed correctly", {
  expect_equal(levels(r$SampleID), c("A.1", "A.2", "B.1", "B.2", "C.1"))
})

test_that("Diversity measurement is parsed correctly", {
  expect_equal(r[1,"diversity"], 5) 
})

context('rarefaction_stats')

s <- rarefaction_stats(r)

test_that("Sample IDs are present in correct column", {
  expect_equal(
    levels(s$SampleID), 
    c("A.1", "A.2", "B.1", "B.2", "C.1"))
})

test_that("Rarefaction depth is present in correct column", {
  expect_equal(unique(s$sequences_per_sample), c(10, 207, 404))
})

test_that("Summary columns are present with correct names", {
  expect_equal(length(s$diversity.mean), 15)
})
