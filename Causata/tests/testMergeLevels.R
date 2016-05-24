# Test of the MergeLevels class
# 
# Author: jasonm
###############################################################################
library(testthat)
library(Causata)

context("MergeLevels")

numbers <- factor(c(
  "one",
  "two","two",
  "three","three","three",
  "four","four","four","four",
  "five","five","five","five","five",
  "six","six","six","six","six","six"))

num.original.factors <- length(levels(numbers))

#test that you get the number of levels your request
#test that you get an 'other' category where appropriate

#iterate from 2 levels all the way to the original number of levels
for (max.levels in seq(2,num.original.factors)){
  mergedNumbers <- MergeLevels(numbers,max.levels)
  
  #test that we have the expected number of levels after merging (including if we ask for the same number we start with
  test_that(paste("Correctly reduce to",max.levels,"levels"), 
            expect_equal(max.levels, length(levels(mergedNumbers))) 
           )
  
  #test that if we have reduced the number of levels, an 'Other' category has appeared, but not otherwise
  if (max.levels < num.original.factors){
      test_that("Other category has been created when categories merged", 
                expect_that("Other" %in% levels(mergedNumbers), is_true())          
      )
  }
  else {
    test_that("Other category has not been created", 
              expect_that("Other" %in% levels(mergedNumbers), is_false())          
    )
  }
}

#Check that you get an error if requesting fewer than 2 levels
expect_that(MergeLevels(numbers,0), throws_error())
expect_that(MergeLevels(numbers,1), throws_error())

#Check that you get an error if you call MergeLevels on a variable that is not a factor
expect_that(MergeLevels(c("a","b","c"),2), throws_error())
