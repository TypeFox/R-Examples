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
  names(sdat_labels) <- names(sdat)
  attributes(sdat)$variable.labels <- sdat_labels
}

#s <- as.surveydata(sdat)
#is.surveydata(s)
#varlabels(s)
#identical(varlabels(s), sdat_labels)

#------------------------------------------------------------------------------

context("Surveydata")

test_that("as.surveydata and is.surveydata works as expected", {
      s <- as.surveydata(sdat)
      #expected_pattern <- c("^", "(_[[:digit:]])*(_.*)?$")
      expected_pattern <- list(sep="_", exclude="other")
      expect_that(s, is_a("surveydata"))
      expect_that(s, is_a("data.frame"))
      expect_that(is.surveydata(s), is_true())
      expect_that(is.surveydata(sdat), is_false())
      expect_that(pattern(s), equals(expected_pattern))
      
      #new_pattern <- c("", "new_pattern$")
      new_pattern <- list(sep=":", exclude="last") 
      s <- as.surveydata(sdat, ptn=new_pattern)
      expect_that(s, is_a("surveydata"))
      expect_that(is.surveydata(s), is_true())
      expect_that(pattern(s), equals(new_pattern))
    })

test_that("Varlabel names are allocated correctly",{
      tdat <- sdat
      attributes(tdat)$variable.labels <- unname(attributes(tdat)$variable.labels)
      t <- as.surveydata(sdat)
      expect_equal(names(t), names(varlabels(t)))
    })

#------------------------------------------------------------------------------

test_that("Varlabel functions work as expected", {
      s <- as.surveydata(sdat)
      expect_equal(varlabels(s), sdat_labels)
      
      varlabels(s) <- 1:8
      expect_equal(varlabels(s), 1:8)
      
      varlabels(s)[3] <- 20
      expect_equal(varlabels(s), c(1:2, 20, 4:8))
      
      s <- as.surveydata(sdat)
      varlabels(s)["crossbreak"] <- "New crossbreak"
      retn <- sdat_labels
      retn["crossbreak"] <- "New crossbreak"
      expect_equal(varlabels(s), retn)
      
    })

#------------------------------------------------------------------------------

test_that("Pattern functions work as expected", {
      pattern <- "-pattern-"
      s <- as.surveydata(sdat)
      attr(s, "pattern") <- pattern
      expect_that(pattern(s), equals(pattern))
      
      attr(s, "pattern") <- NULL
      expect_that(is.null(pattern(s)), is_true())
      pattern(s) <- pattern
      expect_that(attr(s, "pattern"), equals(pattern))
    })

test_that("Removing attributes work as expected", {
      s <- as.surveydata(sdat)
      
      t <- rm.attrs(s)      
      expect_equal(varlabels(t), NULL)
      expect_equal(pattern(t), NULL)
      
      t <- as.data.frame(s, rm.pattern=TRUE)
      expect_equal(t, sdat)

    })

#------------------------------------------------------------------------------

test_that("Name_replace works as expected", {
      s <- as.surveydata(sdat)
      spat <- pattern(s)
      
      names(s) <- gsub("id", "RespID", names(s))
      expect_equal(names(s)[1], "RespID")
      expect_equal(names(varlabels(s))[1], "RespID")
      expect_equal(pattern(s), spat)
      
      newpat <- c("X", "Y")
      s <- as.surveydata(sdat, ptn=newpat)
      
      names(s) <- gsub("id", "RespID", names(s))
      expect_equal(names(s), c("RespID", names(s)[-1]))
      expect_equal(names(varlabels(s)), c("RespID", names(s)[-1]))
      expect_equal(pattern(s), newpat)
      
    })

#------------------------------------------------------------------------------

test_that("warnings are issued when names and varlabels mismatch", {
      s <- as.surveydata(sdat)
      
      s2 <- s
      varlabels(s2) <- varlabels(s2)[-1]
      shows_message(is.surveydata(s2))
      
    })



