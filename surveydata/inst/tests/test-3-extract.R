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

context("Extract")

test_that("`$` extracts correct columns", {
      s <- as.surveydata(sdat, renameVarlabels=TRUE)
      
      expect_equal(s$id, 1:4)
      expect_equal(s$Q4_1, c(1, 2, 1, 2))
      expect_is(s, "surveydata")
    })




#------------------------------------------------------------------------------

#context("Surveydata `[` simple extract")

test_that("`[` simple extract returns surveydata object", {
      s <- as.surveydata(sdat, renameVarlabels=TRUE)
#      load_all(pkg)
#      which.q(s, "Q4")
#      x <- NULL
#      x <- s[2]
#      x
#      varlabels(x)
#      str(x)
      
      expect_is(s[], "surveydata")
      expect_is(s[, 2], "surveydata")
      expect_is(s[1, ], "surveydata")
      expect_is(s[2, 2], "surveydata")
      expect_is(s[, "Q1"], "surveydata")
      expect_is(s[, "Q4"], "surveydata")
    })

test_that("`[` simple extract returns correct data", {
      s <- as.surveydata(sdat, renameVarlabels=TRUE)
      
      expect_equal(s[], s)
      expect_equal(rm.ca(s[2, ]), rm.ca(sdat[2, ]))
      expect_equal(rm.ca(s[, 2]), rm.ca(sdat[, 2, drop=FALSE]))
      expect_equal(rm.ca(s[, "Q1"]), rm.ca(sdat[, 2, drop=FALSE]))
      expect_equal(rm.ca(s[, "Q4"]), rm.ca(sdat[, 3:5, drop=FALSE]))
      expect_equal(rm.ca(s[2, "Q4"]), rm.ca(sdat[2, 3:5, drop=FALSE]))
      expect_equal(rm.ca(s[1:2, "Q10"]), rm.ca(sdat[1:2, 7, drop=FALSE]))
      expect_equal(rm.ca(s[, "weight"]), rm.ca(sdat[, "weight", drop=FALSE]))
    })

test_that("`[` simple extract returns correct varlabels", {
      s <- as.surveydata(sdat, renameVarlabels=TRUE)
      
      expect_equal(varlabels(s[]), sdat_labels)
      expect_equal(varlabels(s[2]), sdat_labels[2])
      expect_equal(varlabels(s[, 2]), sdat_labels[2])
      expect_equal(varlabels(s[2:4, 5]), sdat_labels[5])
      expect_equal(varlabels(s[, "Q1"]), sdat_labels[2])
      expect_equal(varlabels(s[, "Q4"]), sdat_labels[3:5])
      expect_equal(varlabels(s[2, "Q4"]), sdat_labels[3:5])
      expect_equal(varlabels(s[2, "Q1"]), sdat_labels[2])
      
    })

test_that("`[` simple extract with drop=TRUE returns vectors", {
  s <- as.surveydata(sdat, renameVarlabels=TRUE)
  
  expect_equal(rm.ca(s[, 2, drop=TRUE]), rm.ca(sdat[, 2, drop=TRUE]))
  expect_equal(rm.ca(s[, "Q1", drop=TRUE]), rm.ca(sdat[, 2, drop=TRUE]))
  expect_equal(rm.ca(s[, "Q4", drop=TRUE]), rm.ca(sdat[, 3:5, drop=TRUE]))
  expect_equal(rm.ca(s[2, "Q4", drop=TRUE]), rm.ca(sdat[2, 3:5, drop=TRUE]))
  expect_equal(rm.ca(s[1:2, "Q10", drop=TRUE]), rm.ca(sdat[1:2, 7, drop=TRUE]))
  expect_equal(rm.ca(s[, "weight", drop=TRUE]), rm.ca(sdat[, "weight", drop=TRUE]))
  
})


#------------------------------------------------------------------------------
    
#context("Surveydata `[` complex extract")
test_that("`[` complex extract works as expected", {
      s <- as.surveydata(sdat, renameVarlabels=TRUE)
      
      expect_equal(rm.ca(s[, c(1, 3)]), rm.ca(sdat[, c(1, 3)]))
      expect_equal(rm.ca(s[, -1]), rm.ca(sdat[, -1]))
      expect_equal(rm.ca(s[, c(1, "Q4")]), rm.ca(sdat[, c(1, 3:5)]))
      expect_equal(rm.ca(s[, c("Q1", "Q4")]), rm.ca(sdat[, c(2, 3:5)]))

      expect_equal(varlabels(s[, c(1, 3)]), sdat_labels[c(1, 3)])
      expect_equal(varlabels(s[, -1]), sdat_labels[-1])
      expect_equal(varlabels(s[, c(1, "Q4")]), sdat_labels[c(1, 3:5)])
      expect_equal(varlabels(s[, c("Q1", "Q4")]), sdat_labels[c(2, 3:5)])
      
    })
          

#------------------------------------------------------------------------------

test_that("`[` extract with logicals", {
      s <- as.surveydata(sdat, renameVarlabels=TRUE)
      
      i <- sdat$id==1
      j <- grepl("Q4", names(s))
      expect_equal(rm.ca(s[i, ]), rm.ca(sdat[i, ]))
      expect_equal(rm.ca(s[!i, ]), rm.ca(sdat[!i, ]))
      expect_equal(rm.ca(s[i, j]), rm.ca(sdat[i, j]))
      expect_equal(rm.ca(s[i, !j]), rm.ca(sdat[i, !j]))
      expect_equal(rm.ca(s[!i, j]), rm.ca(sdat[!i, j]))
      expect_equal(rm.ca(s[!i, j]), rm.ca(sdat[!i, j]))
      expect_equal(rm.ca(s[!i, !j]), rm.ca(sdat[!i, !j]))
      
    })
