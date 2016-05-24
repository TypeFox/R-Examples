context("read_SEQ2R")

test_that("file is *.SEQ file",{
  expect_error(read_SEQ2R(file = 2),"file has to be a character")

})
