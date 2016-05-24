# 
# Author: andrie
###############################################################################

context("tools")

test_that("dropout calculation is correct", {
      
      rest <- setNames(c(215, 35), c("id", "Q23_other"))
      expect_equal(dropout(membersurvey[-(108:109)]), rest)
    })

