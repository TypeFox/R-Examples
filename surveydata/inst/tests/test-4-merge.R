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
      Q10 = factor(c("Male", "Female", "Female", "Male")),
      crossbreak  = c("A", "A", "B", "B"), 
      crossbreak2 = c("D", "E", "D", "E"),
      weight      = c(0.9, 1.1, 0.8, 1.2)
  )
  
  sdat_labels <- c(
      "RespID",
      "Question 1", 
      "Question 4: red", "Question 4: green", "Question 4: blue", 
      "Question 10",
      "crossbreak",
      "crossbreak2",
      "weight")
  varlabels(sdat) <- sdat_labels
}



context("Merge")

test_that("Merge of surveyordata objects work as expected",{
      sdat2 <- data.frame(
          id   = 5:6,
          Q1   = c("Yes", "No"),
          Q4_1 = c(5, 6), 
          Q4_2 = c(7, 8), 
          Q4_3 = c(9, 10), 
          crossbreak  = c("U", "V"), 
          crossbreak2 = c("X", "Y"),
          weight      = c(0.95, 1.05)
      )
      varlabels(sdat2) <- names(sdat2)
      
      s1 <- as.surveydata(sdat, renameVarlabels=TRUE)
      s2 <- as.surveydata(sdat2, ptn=c("_", ""), renameVarlabels=TRUE)
      expect_that(
          sm <- merge(s1, s2, all=TRUE),
          gives_warning("In merge of surveydata objects, patterns of objects differ")
      )
      #sm <- merge(s1, s2, all=TRUE)
      expect_that(sm, is_a("surveydata"))
      expect_equal(nrow(sm), 6)
      expect_equal(pattern(s1), pattern(sm))
    })