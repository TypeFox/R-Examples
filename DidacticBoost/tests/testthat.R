Sys.setenv("R_TESTS" = "")
if(require(testthat) & require(rpart)){
  test_check("DidacticBoost")
}
