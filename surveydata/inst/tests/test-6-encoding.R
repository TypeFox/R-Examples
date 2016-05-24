# Unit tests using package testthat
# 
# Author: Andrie
#------------------------------------------------------------------------------

context("Encoding")

test_that("encoding functions work as expected", {

      expect_equal(encToInt("\\xfa"), c(92, 120, 102, 97))
      #expect_equal(intToEnc(8212), "—")
      expect_equal(intToEnc(encToInt("\\xfa")), "\\xfa")
      expect_equal(encToInt(intToEnc(250)), 250)
      
#      print(encToInt(intToEnc(8212)))
#      print(intToEnc(encToInt("\\xfa")))
    })


