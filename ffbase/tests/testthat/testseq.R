library(testthat)
library(ff)

context("ffseq")

oldffmaxbytes <- getOption("ffmaxbytes")
options(ffmaxbytes = 4)

test_that("ffseq_len works",{
	test.ff <- ffseq_len(1000000)
	test.ram <- seq_len(1000000)
	expect_equal(test.ram, test.ff[])
})

test_that("ffseq works",{
  test.ff <- ffseq(0, 1, length.out=1000)
  test.ram <- seq(0, 1, length.out=1000)
  expect_equal(test.ram, test.ff[])
})
test_that("ffseq works",{
  test.ff <- ffseq(stats::rnorm(20))
  test.ram <- seq(stats::rnorm(20))
  expect_equal(test.ram, test.ff[])
})
test_that("ffseq works",{
  test.ff <- ffseq(1, 9, by = 2) 
  test.ram <- seq(1, 9, by = 2) 
  expect_equal(test.ram, test.ff[])
})
test_that("ffseq works",{
  test.ff <- ffseq(1, 9, by = pi)
  test.ram <- seq(1, 9, by = pi) 
  expect_equal(test.ram, test.ff[])
})
test_that("ffseq works",{
  test.ff <- ffseq(1, 6, by = 3)
  test.ram <- seq(1, 6, by = 3)
  expect_equal(test.ram, test.ff[])
})
test_that("ffseq works",{
  test.ff <- ffseq(1.575, 5.125, by=0.05)
  test.ram <- seq(1.575, 5.125, by=0.05)
  expect_equal(test.ram, test.ff[])
})
test_that("ffseq works",{
  test.ff <- ffseq(17) 
  test.ram <- seq(17) 
  expect_equal(test.ram, test.ff[])
})

options(ffmaxbytes = oldffmaxbytes)



