library(testthat)

context("Linear editrow derivarions")

test_that("Various rows work",{
   edt <- parse(text=c("x==y", "x+w ==y"))
   e <- parseNum(edt[[1]])
   #print(e)
})

test_that("Parsing a constant works",{
   edt <- parse(text=c("x < 2"))
   e <- parseNum(edt[[1]])
   #print(e)
})

test_that("Parsing a inequality works",{
   edt <- parse(text=c("x > 2"))
   e <- parseNum(edt[[1]])
   #print(e)
})

test_that("Parsing a negative coefficient works",{
   edt <- parse(text=c("x == -2"))
   #e <- makeEditRow(edt[[1]])
   #print(e)
})
