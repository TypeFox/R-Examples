# Unit tests for "surveydata" class
# 
# Author: Andrie
#------------------------------------------------------------------------------


# Reads sample data obtained from http://www3.norc.org/GSS+Website/Download/SPSS+Format/

#library(foreign)
#path <- file.path("f:", "git", "surveydata", "surveydata", "inst", "tests")
#filename <- "2010.sav"
#gss <- read.spss(file=file.path(path, filename), to.data.frame=TRUE)
#save(gss, file=file.path(path, "gss.rda"))

path <- file.path("inst", "tests")
filename <- "gss.rda"
load(file=file.path(path, filename))

gss <- as.surveydata(gss, ptn=c("^", "(.*?)$"))

context("Integrated test of reading and analysing a survey file")

  test_that("surveydata works with real spss data file", {
        
    qs <- questions(gss)
    
    expect_equal(length(qs), 790)
    expect_equal(head(qs, 20), 
        c("mar1", "mar2", "mar3", "mar4", "mar5", "mar6", "mar7", "mar8", 
        "mar9", "mar11", "mar12", "abany", "abdefect", "abhlth", "abnomore", 
        "abpoor", "abrape", "absingle", "acqntsex", "adults"))
    
    expect_equal(head(varlabels(gss)), 
      structure(c("MARITAL STATUS OF 1ST PERSON", "MARITAL STATUS OF 2ND PERSON", 
              "MARITAL STATUS OF 3RD PERSON", "MARITAL STATUS OF 4TH PERSON", 
              "MARITAL STATUS OF 5TH PERSON", "MARITAL STATUS OF 6TH PERSON"
          ), .Names = c("mar1", "mar2", "mar3", "mar4", "mar5", "mar6")))
  
    expect_equal((qText(gss, "where")), 
      c("WHERE IS 1ST PERSON STAYING?", "WHERE IS 11TH PERSON (VISITOR) STAYING?", 
          "WHERE IS 2ND PERSON STAYING?", "WHERE IS 3RD PERSON STAYING?", 
          "WHERE IS 4TH PERSON STAYING?", "WHERE IS 5TH PERSON STAYING?", 
          "WHERE IS 6TH PERSON STAYING?", "WHERE IS 7TH PERSON STAYING?"))

})

