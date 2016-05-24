context('read_qiime_distmat')

d <- read_qiime_distmat(
  system.file("testdata", "distmat.txt", package="qiimer"))

test_that("A dist object is produced", {
  expect_is(d, "dist")
})

test_that("Labels are correct", {
  expect_equal(attr(d, "Labels"), c(
    "A1002.YHNS", "A1005.HELV", "A1006.KQJA", "A1012.AOPN", "A1013.FFXK", 
    "A1015.FREJ", "A1016.BWDG", "A1038.SYHH", "A1043.BVIY", "H2O.1", 
    "SG.H2O.2"))
})

test_that("Values are correct", {
  expect_equal(d[1], 0.8749526)
})

dm_names <- letters[1:4]
dm <- matrix(
  c(0, 1, 2, 3,
    1, 0, 4, 5,
    2, 4, 0, 6,
    3, 5, 6, 0), 
  ncol=4, 
  dimnames=list(dm_names, dm_names))

context('dist_get')

test_that('dist_get works with named indices', {
  expect_equal(
    dist_get(dm, c("a", "b"), c("c", "d")), c(2, 5))
  expect_equal(
    dist_get(dm, c("a", "a", "a"), c("a", "b", "c")), c(0, 1, 2))
})

test_that('dist_get recycles short vectors', {
  expect_equal(dist_get(dm, "a", c("a", "b", "c")), c(0, 1, 2))
})

test_that('dist_get works with numeric indices', {
  expect_equal(dist_get(dm, c(1, 3), c(4, 4)), c(3, 6))
})

context('dist_subset')

test_that('dist_subset returns a dist object', {
  expect_equal(class(dist_subset(dm, 1:3)), "dist")
})

context("dist_groups")

test_that("dist_groups labels groups correctly", {
  dg <- dist_groups(dm, c("A", "A", "B", "B"))
  expect_equal(levels(dg$Label), c("Between A and B", "Within A", "Within B"))
  dg <- dist_groups(dm, c("A", "B", "A", "B"))
  expect_equal(levels(dg$Label), c("Between A and B", "Within A", "Within B"))
})
