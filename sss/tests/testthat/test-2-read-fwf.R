# Add comment
# 
# Author: Andrie
#----------------------------------------------------------------------------------

if(interactive()){
  library(testthat)
  sampleRoot <- "sss/tests/testthat/samples/sample-0"
} else {
  sampleRoot <- "samples/sample-0"
}

context("Read fwf")


test_that("widths work", {
      ff <- tempfile()
      cat(file = ff, "123456", "987654", sep="\n")
      
      rest <- structure(list(
              V1 = c(1, 9),
              V2 = c(23, 87),
              V3 = c(456, 654)),
            .Names = c("V1", "V2", "V3"),
            row.names = c(NA, -2L),
            class="data.frame"
            )
      
      test <- fast.read.fwf(ff, widths = c(1, 2, 3))    #> 1 23 456 \ 9 87 654
      expect_equal(rest, test)
      expect_is(test, "data.frame")
      unlink(ff)
    
    })

test_that("colclasses work",{
      ff <- tempfile()
      cat(file=ff, "abc123def1010", "ghi456jkl1100", sep="\n")
      
      test <- fast.read.fwf(file=ff, widths=c(3, 3, 3, 1, 1, 1, 1),
          colClasses=c("character", "numeric", "character", rep("logical", 4)))
      rest <- structure(list(
              V1 = c("abc", "ghi"), 
              V2 = c(123, 456), 
              V3 = c("def", "jkl"), 
              V4 = c(TRUE, TRUE), 
              V5 = c(FALSE, TRUE), 
              V6 = c(TRUE, FALSE), 
              V7 = c(FALSE, FALSE)), 
          .Names = c("V1", "V2", "V3", "V4", "V5", "V6", "V7"), 
          row.names = c(NA, -2L), class = "data.frame")
      unlink(ff)
      expect_equal(test, rest)
    })



