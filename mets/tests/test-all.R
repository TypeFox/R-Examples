if (require(testthat)) {
  library(mets)
  if (exists("test_check"))
      test_check("mets")
}
