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
      Q4_other = LETTERS[1:4],
      Q10 = factor(c("Male", "Female", "Female", "Male")),
      crossbreak  = c("A", "A", "B", "B"), 
      crossbreak2 = c("D", "E", "D", "E"),
      weight      = c(0.9, 1.1, 0.8, 1.2)
  )
  
  sdat_labels <- c(
      "RespID",
      "Question 1", 
      "Question 4: red", "Question 4: green", "Question 4: blue", "Question 4: other",
      "Question 10",
      "crossbreak",
      "crossbreak2",
      "weight"
  )
  names(sdat_labels) <- names(sdat)
  varlabels(sdat) <- sdat_labels
}

rm.ca <- function(x){
  class(x) <- class(x)[!grepl("surveydata", class(x))]
  rm.attrs(x)
}

#------------------------------------------------------------------------------

context("Replace")

test_that("`$<-` NULL removes column as well as varlabel", {
      s <- as.surveydata(sdat, renameVarlabels=TRUE)
      s$id <- NULL
      expect_true(is.na(match("id", names(s))))
      expect_true(is.na(match("id", names(varlabels(s)))))
      expect_equal(names(s), names(sdat[-1]))
      expect_equal(names(varlabels(s)), names(sdat[-1]))
      expect_is(s, "surveydata")
    })

test_that("`$<-` existing_name maintains correct varlabels",{
      s <- as.surveydata(sdat, renameVarlabels=TRUE)
      expect_equal(varlabels(sdat), varlabels(s))
      s$Q4_1 <- 1:4
      expect_equal(varlabels(sdat), varlabels(s))
    })

test_that("`$<-` newname inserts column and new varlabel", {
      s <- as.surveydata(sdat, renameVarlabels=TRUE)
      s$newid <- 101:104
      expect_equal(s$newid, 101:104)
      expect_true(all(s$newid==101:104))
      expect_false(is.na(match("newid", names(varlabels(s)))))
      expect_is(s, "surveydata")
    })

#------------------------------------------------------------------------------

test_that("`[<-` NULL removes column as well as varlabel", {
      s <- as.surveydata(sdat, renameVarlabels=TRUE)
      s[, "id"] <- NULL
      #browser()
      expect_true(is.na(match("id", names(s))))
      expect_true(is.na(match("id", names(varlabels(s)))))
      expect_equal(names(s), names(sdat[-1]))
      expect_equal(names(varlabels(s)), names(sdat[-1]))
      expect_is(s, "surveydata")
    })

test_that("`[<-` existing_name maintains correct varlabels",{
      s <- as.surveydata(sdat, renameVarlabels=TRUE)
      expect_equal(varlabels(sdat), varlabels(s))
      s[, "Q4_1"] <- 1:4
      expect_equal(varlabels(sdat), varlabels(s))
    })

#------------------------------------------------------------------------------

test_that("`[<-` newname inserts column and new varlabel", {
      s <- as.surveydata(sdat, renameVarlabels=TRUE)
      s["newid"] <- 101:104
      expect_equal(s$newid, 101:104)
      expect_false(is.na(match("newid", names(varlabels(s)))))
      expect_is(s, "surveydata")
    })

