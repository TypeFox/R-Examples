

context("adjustRecords")

test_that("weights are used",{
  E <- editmatrix(expression(
    x + y == 1,
    x > 0,
    y >= 0)
  )
  dat <- data.frame(x=1,y=2)
  expect_equal(adjustRecords(E,dat)$adjusted, data.frame(x=0,y=1),tolerance=0.01)
  expect_equal(
    unlist(adjustRecords(E,dat,w=c(2,1))$adjusted)
    , adjust(E,unlist(dat[1,]),w=c(2,1))$x 
  )
  # extra, unrelated variable
  dat <- data.frame(x=1,y=2,z=0)
  expect_equal(adjustRecords(E,dat)$adjusted, data.frame(x=0,y=1,z=0),tolerance=0.01)
  expect_equal(
    unlist(adjustRecords(E,dat,w=c(2,1,1))$adjusted)
    , adjust(E,unlist(dat[1,]),w=c(2,1,1))$x 
  )
})


