
library(testthat)

context("Parcels")

test_that("wrapping and contents work", {
  abc <- letters[1:3]
  pcl <- wrap(abc)
  expect_equivalent(contents(pcl), letters[1:3])
  abc <- letters[4:6]
  expect_equivalent(contents(pcl), letters[4:6])
  
  names(abc) <- LETTERS[1:3]
  expect_equivalent(names(contents(pcl)), LETTERS[1:3])
})

test_that("wrapping into an environment works", {
  env <- new.env(emptyenv())
  env$abc <- letters[1:3]
  pcl <- wrap(abc, env=env)
  expect_equivalent(contents(pcl), letters[1:3])
  env$abc <- letters[4:6]
  expect_equivalent(contents(pcl), letters[4:6])
})

test_that("wrapping arbitrary expressions works", {
  pcl <- wrap(x^3)
  x <- 2
  expect_equivalent(contents(pcl), 8)
  x <- 3
  expect_equivalent(contents(pcl), 27)
})

test_that("unwrap_as works", {
  abc <- letters[1:3]
  pcl <- wrap(abc)
  unwrap_as(ref, pcl)
  expect_equivalent(ref, abc)
  abc <- 1:3
  expect_equivalent(ref, abc)
})

test_that("assigning to contents works", {
  abc <- letters[1:3]
  pcl <- wrapset(abc,)
  contents(pcl) <- letters[24:26]
  expect_equivalent(abc, letters[24:26])
  expect_equivalent(contents(pcl), letters[24:26])
})

