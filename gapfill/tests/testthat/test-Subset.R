#require(testthat); source("../../gapfill.R")


context("test-Subset")

test_that("ArrayAround-trivial",{
  data <- array(1, c(1,1,1,1))
  out <- ArrayAround(data = data,
                     mp = c(1L,1L,1L,1L),
                     size = c(0L,0L,0L,0L))
  expect_identical(attr(out, "mp"), c(1L, 1L, 1L, 1L))
  attr(out, "mp") <- NULL
  expect_identical(out, data)

  expect_error(ArrayAround(data = array(1:3, c(1,1,1)),
                           mp = c(1,1,1,1),
                           size = c(0,0,0,0)))
  expect_error(ArrayAround(data = data,
                           mp = c(1,1,1,1),
                           size = c(-1,0,0,0)))
  expect_error(ArrayAround(data = data,
                           mp = c(1,1,1,-1),
                           size = c(1,0,0,0)))
  expect_error(ArrayAround(data = data,
                           mp = c(1,1,1,2),
                           size = c(1,0,0,0)))

})

test_that("ArrayAround-non-trivial",{
  data <- array(1:24, 1:4)
  for(i in 1:24){
      i4 <- IndexOneFour(10, dim(data))
      expect_equal(c(data[i4[1],i4[2],i4[3],i4[4]]),
                   c(ArrayAround(data, i4, c(0,0,0,0))))
      expect_equal(data[1,i4[2],i4[3],i4[4]],
                   c(ArrayAround(data, i4, c(1,0,0,0))))
      expect_equal(c(data[1,1:2,i4[3],i4[4]]),
                   c(ArrayAround(data, i4, c(1,2,0,0))))
      expect_equal(c(data[1,1:2,1:3,i4[4]]),
                   c(ArrayAround(data, i4, c(1,2,3,0))))
      expect_equal(c(data[1,1:2,1:3,1:4]),
                   c(ArrayAround(data, i4, c(1,2,3,4))))
  }
})


test_that("ArrayAroundRandom",{
  expect_error(ArrayAroundRandom(data = array(1:24, c(2,3,3,2)),
                                 target = "missing"))
  expect_error(ArrayAroundRandom(data = array(NA, c(2,3,3,2)),
                                 target = "observed"))
  out <- ArrayAroundRandom(data = array(NA, c(2,3,3,2)),
                           target = "all")
  expect_false(is.null(attr(out, "mp")))
  expect_true(is.array(out))
  expect_equal(dim(out), c(1,1,1,1))
})


test_that("Subset",{
  data <- array(1:24, c(2,4,3,2))
  for(i in 0:3){
      expect_equal(c(Subset(data, mp = c(1, 1, 1, 1), i = i,
                            initialSize = c(0, 0, 0, 0))),
                   c(data[1:min(2,i+1),1:min(4,i+1),1,1]))
  }
})



