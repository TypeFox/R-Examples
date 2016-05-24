# Unit tests for "surveydata" class
# 
# Author: Andrie
#------------------------------------------------------------------------------

{
  sdat <- data.frame(
      id   = 1:4,
      Q1   = c("Yes", "No", "Yes", "Yes"),
      Q4_1 = c(1, 2, 1, 2), 
      Q4_2 = c(3, 4, 4, 3), 
      Q4_3 = c(5, 5, 6, 6), 
      Q4_other = c(NA, NA, "some text", NA),
      Q10 = factor(c("Male", "Female", "Female", "Male")),
      crossbreak  = c("A", "A", "B", "B"), 
      crossbreak2 = c("D", "E", "D", "E"),
      weight      = c(0.9, 1.1, 0.8, 1.2)
  )
  
  varlabels(sdat) <- c(
      "RespID",
      "Question 1", 
      "Question 4: red", "Question 4: green", "Question 4: blue", 
      "Question 4: Other",
      "Question 10",
      "crossbreak",
      "crossbreak2",
      "weight")
  
  sdat2 <- data.frame(
      id   = 1:4,
      Q1   = c("Yes", "No", "Yes", "Yes"),
      `Q4__1` = c(1, 2, 1, 2), 
      `Q4__2` = c(3, 4, 4, 3), 
      `Q4__3` = c(5, 5, 6, 6), 
      `Q4__ignore` = c(NA, NA, "some text", NA),
      Q10 = factor(c("Male", "Female", "Female", "Male")),
      crossbreak  = c("A", "A", "B", "B"), 
      crossbreak2 = c("D", "E", "D", "E"),
      weight      = c(0.9, 1.1, 0.8, 1.2),
      check.names=FALSE
  )
  
  varlabels(sdat2) <- c(
      "RespID",
      "Question 1", 
      "Question 4: red", "Question 4: green", "Question 4: blue", 
      "Question 4: Other",
      "Question 10",
      "crossbreak",
      "crossbreak2",
      "weight")
  
}

rm.ca <- function(x){
  class(x) <- class(x)[!grepl("surveydata", class(x))]
  rm.attrs(x)
}

#------------------------------------------------------------------------------

context("which.q")

test_that("which.q returns correct question positions", {
      s <- as.surveydata(sdat, renameVarlabels=TRUE)
      expect_that(which.q(s, c(1)), equals(1))
      expect_that(which.q(s, c(4)), equals(4))
      expect_that(which.q(s, c(-1)), equals(-1))
      expect_that(which.q(s, "Q1"), equals(2))
      expect_that(which.q(s, "Q10"), equals(7))
      expect_that(which.q(s, "Q4"), equals(3:5))
      expect_that(which.q(s, "Q2"), equals(integer(0)))
      
      expect_that(which.q(s, c("Q1", "Q4")), equals(c(2, 3:5)))
      expect_that(which.q(s, c("Q1", "crossbreak")), equals(c(2, 8)))
      expect_that(which.q(s, c("Q4", "crossbreak2")), equals(c(3:5, 9)))
      
      expect_that(which.q(s, c(3, "crossbreak2")), equals(c(3, 9)))
      
    })

#------------------------------------------------------------------------------

#context("which.q 2")

test_that("which.q returns correct question positions", {
      s2 <- as.surveydata(sdat2, ptn=list(sep="__", exclude="ignore"), renameVarlabels=TRUE)
      expect_that(which.q(s2, c(1)), equals(1))
      expect_that(which.q(s2, c(4)), equals(4))
      expect_that(which.q(s2, c(-1)), equals(-1))
      expect_that(which.q(s2, "Q1"), equals(2))
      expect_that(which.q(s2, "Q10"), equals(7))
      expect_that(which.q(s2, "Q4"), equals(3:5))
      expect_that(which.q(s2, "Q2"), equals(integer(0)))
      
      expect_that(which.q(s2, c("Q1", "Q4")), equals(c(2, 3:5)))
      expect_that(which.q(s2, c("Q1", "crossbreak")), equals(c(2, 8)))
      expect_that(which.q(s2, c("Q4", "crossbreak2")), equals(c(3:5, 9)))
      
      expect_that(which.q(s2, c(3, "crossbreak2")), equals(c(3, 9)))
      
    })


