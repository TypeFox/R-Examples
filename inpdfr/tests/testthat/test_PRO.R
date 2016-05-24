library(testthat)
library(inpdfr)

context("getListFiles")

test_that("getListFiles returns a list of length 2",{
  expect_is(getListFiles(mywd=getwd()), "list")
  expect_equal(length(getListFiles(mywd=getwd())), 2)
})

context("quitSpaceFromChars")

test_that("quitSpaceFromChars works with all types of file names",{
  expect_is(quitSpaceFromChars(vectxt="file.pdf"),"logical")
  expect_is(quitSpaceFromChars(vectxt="f i l e.pdf"),"logical")
  expect_is(quitSpaceFromChars(vectxt="f,i.l-e.pdf"),"logical")
  expect_is(quitSpaceFromChars(vectxt="f%i'le.pdf"),"logical")
})

